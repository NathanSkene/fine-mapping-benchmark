"""Eager-reference equivalence for complete BEATRICE candidate rounds."""
import copy
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import random
import sys
import time
import pickle
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "batch_posterior_memo", ROOT / "BEATRICE_annot_sparse/scripts/posterior_memo.py")
pm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pm)


def record_measurement(value):
    print("BATCH_EQUIVALENCE_MEASUREMENT " + json.dumps(value, sort_keys=True))
    path = os.getenv("FMB_BATCH_METRICS_PATH")
    if path:
        path = Path(path)
        records = json.loads(path.read_text()) if path.exists() else []
        records.append(value)
        path.write_text(json.dumps(records, indent=2) + "\n")


class ConfigurationBatchTests(unittest.TestCase):
    def context(self, prior=(.2, .4, .8)):
        memo = {(0,): 0.}
        context = pm.freeze_memo(memo, prior)
        return memo, context

    def test_eager_reference_exact_support_order_weights_and_null(self):
        for p in (3, 12, 300):
            rng = random.Random(19 + p)
            prior = [rng.uniform(.01, .9) for _ in range(p)]
            updates = [((j,), rng.uniform(-10, 10)) for j in range(p)]
            updates += [(tuple(sorted(rng.sample(range(p), 2))), rng.uniform(-10, 10))
                        for _ in range(p)]
            updates += [(updates[0][0], 0.), ((), 0.), (updates[1][0], -10000.)]
            eager, original = self.context(prior)
            batched, context = copy.deepcopy((eager, original))
            for key, log_bf in updates:
                pm.add_configuration(eager, original, key, log_bf)
            with pm.configuration_batch(batched, context) as batch:
                for key, log_bf in updates:
                    self.assertIsNone(batch.add(key, log_bf))
                    self.assertTrue(batch.contains(key))
            self.assertEqual(context, original)
            self.assertEqual(list(context["log_bf"]), list(original["log_bf"]))
            self.assertEqual(list(batched), list(eager))
            self.assertEqual(batched, eager)
            self.assertEqual(batched[()], eager[()])
            self.assertEqual(context["final_prior"], prior)

    def test_boundary_prior_and_zero_log_bf_preserve_impossible_models(self):
        prior = [0., 1., .2]
        eager = {(1,): 0.}
        original = pm.freeze_memo(eager, prior)
        batched, context = copy.deepcopy((eager, original))
        updates = [(key, 0.) for k in range(4) for key in itertools.combinations(range(3), k)]
        for key, log_bf in updates:
            pm.add_configuration(eager, original, key, log_bf)
        with pm.configuration_batch(batched, context) as batch:
            for key, log_bf in updates:
                batch.add(key, log_bf)
        self.assertEqual(batched, eager)
        self.assertEqual(batched[()], 0.)
        self.assertEqual(batched[(0, 1)], 0.)
        self.assertAlmostEqual(batched[(1,)], .8)
        self.assertAlmostEqual(batched[(1, 2)], .2)

    def test_bad_additions_reject_immediately_without_publishing_partial_state(self):
        bad = [((0,), float("nan")), ((0,), float("inf")), ((0,), -float("inf")),
               ((), 1.), ((0, 0), 0.), ((1, 0), 0.), ((3,), 0.), ((-1,), 0.), ((True,), 0.)]
        for key, log_bf in bad:
            with self.subTest(key=key, log_bf=log_bf):
                memo, context = self.context()
                initial = copy.deepcopy((memo, context))
                with self.assertRaises(ValueError):
                    with pm.configuration_batch(memo, context) as batch:
                        batch.add((2,), 2.)
                        batch.add(key, log_bf)
                self.assertEqual((memo, context), initial)

    def test_corrupt_context_is_rejected_before_any_round(self):
        for field, value in [("format", "legacy"), ("log_bf", {(0,): float("nan")}),
                             ("log_bf", {(3,): 0.}), ("log_bf", {(): 1.}),
                             ("final_prior", [float("nan"), .4, .8])]:
            memo, context = self.context()
            context[field] = value
            initial_memo = dict(memo)
            with self.assertRaises(ValueError):
                with pm.configuration_batch(memo, context):
                    self.fail("Corrupt context entered a batch")
            self.assertEqual(memo, initial_memo)

    def test_frozen_prior_cache_and_public_memo_cannot_change_during_round(self):
        for mutation in [lambda m, c: c["final_prior"].__setitem__(0, .3),
                         lambda m, c: c["log_bf"].__setitem__((0,), 3.),
                         lambda m, c: m.__setitem__((0,), .999)]:
            memo, context = self.context()
            with self.assertRaisesRegex(ValueError, "changed during"):
                with pm.configuration_batch(memo, context) as batch:
                    batch.add((2,), 2.)
                    mutation(memo, context)
            self.assertNotIn((2,), context["log_bf"])

    def test_empty_aborted_and_reused_batches_have_explicit_boundaries(self):
        memo, context = self.context()
        initial = copy.deepcopy((memo, context))
        with pm.configuration_batch(memo, context) as batch:
            self.assertTrue(batch.contains(()))
        self.assertEqual((memo, context), initial)
        with self.assertRaisesRegex(RuntimeError, "active batch"):
            batch.add((2,), 0.)
        with self.assertRaisesRegex(RuntimeError, "cannot be reused"):
            with batch:
                pass
        with self.assertRaisesRegex(RuntimeError, "interrupt"):
            with pm.configuration_batch(memo, context) as batch:
                batch.add((2,), 0.)
                raise RuntimeError("interrupt")
        self.assertEqual((memo, context), initial)

    def test_moderate_round_normalizes_once_instead_of_per_addition(self):
        p = 500
        prior = [.01 + .02 * j / p for j in range(p)]
        updates = [((j,), math.sin(j / 7)) for j in range(p)]
        eager, original = self.context(prior)
        batched, context = copy.deepcopy((eager, original))
        with mock.patch.object(pm, "posterior_weights", wraps=pm.posterior_weights) as weights, \
                mock.patch.object(pm, "validate_context", wraps=pm.validate_context) as validation:
            began = time.perf_counter()
            for key, log_bf in updates:
                pm.add_configuration(eager, original, key, log_bf)
            eager_seconds = time.perf_counter() - began
            eager_passes = weights.call_count
            eager_validations = validation.call_count
        with mock.patch.object(pm, "posterior_weights", wraps=pm.posterior_weights) as weights, \
                mock.patch.object(pm, "validate_context", wraps=pm.validate_context) as validation:
            began = time.perf_counter()
            with pm.configuration_batch(batched, context) as batch:
                for key, log_bf in updates:
                    batch.add(key, log_bf)
            batch_seconds = time.perf_counter() - began
            batch_passes = weights.call_count
            batch_validations = validation.call_count
        self.assertEqual(eager_passes, p)
        self.assertEqual(batch_passes, 1)
        self.assertEqual(eager_validations, 2 * p)
        self.assertEqual(batch_validations, 2)
        self.assertEqual(batched, eager)
        record_measurement(dict(scope="memo_only_round", p=p, additions=p,
            eager_normalizations=eager_passes, batch_normalizations=batch_passes,
            eager_context_validations=eager_validations, batch_context_validations=batch_validations,
            eager_seconds=eager_seconds, batch_seconds=batch_seconds,
            support_size=len(batched), exact_weights=True))


@unittest.skipUnless(os.getenv("FMB_NEURAL_SMOKE") == "1", "optional real Torch credible-set expansion")
class CredibleSetBatchEquivalenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(ROOT / "BEATRICE_annot_sparse"))
        import numpy as np
        import torch
        from scripts import generate_credible_sets as gen, posterior_memo as native_pm
        cls.np, cls.torch, cls.gen, cls.native_pm = np, torch, gen, native_pm
        torch.set_num_threads(1)

    def expansion(self, p, seed, eager, allow_dup=False):
        np, torch, gen, native_pm = self.np, self.torch, self.gen, self.native_pm
        rng = np.random.default_rng(seed)
        z = torch.tensor(rng.normal(0, .25, p), dtype=torch.float32)
        z[0], z[p // 2] = 8., 6.
        idx = torch.arange(p)
        ld = (.4 ** torch.abs(idx[:, None] - idx[None, :])).float()
        s = torch.matmul(torch.inverse(ld), z.unsqueeze(1))
        p0 = torch.tensor(rng.uniform(.01, .04, p), dtype=torch.float32)
        memo = {(): 0.}
        context = native_pm.freeze_memo(memo, p0)
        # Seed genuine cached likelihoods with the unchanged single-add API.
        gen.abf(z, ld, memo, context, 100, .05, p0, [0], s)
        gen.abf(z, ld, memo, context, 100, .05, p0, [p // 2], s)
        decisions = []
        originals = {name: getattr(gen, name) for name in
            ("find_cond_prob_constrained_causal", "find_cond_prob_constrained_causal_no_dup")}

        class EagerReferenceBatch:
            def __init__(self, memo, context):
                self.memo, self.context = memo, context
            def __enter__(self):
                return self
            def __exit__(self, *args):
                return False
            def contains(self, key):
                return key in self.memo
            def add(self, key, log_bf):
                native_pm.add_configuration(self.memo, self.context, key, log_bf)

        def traced(name):
            def call(*args):
                result = originals[name](*args)
                decisions.append((name, tuple(args[0]), int(result[0]), float(result[1]),
                                  result[2].copy(), tuple(args[2]), dict(args[2])))
                return result
            return call

        with mock.patch.object(gen, "configuration_batch", EagerReferenceBatch if eager
                               else native_pm.configuration_batch), \
                mock.patch.object(gen, "find_cond_prob_constrained_causal", traced("find_cond_prob_constrained_causal")), \
                mock.patch.object(gen, "find_cond_prob_constrained_causal_no_dup", traced("find_cond_prob_constrained_causal_no_dup")), \
                mock.patch.object(native_pm, "posterior_weights", wraps=native_pm.posterior_weights) as weights:
            began = time.perf_counter()
            start = gen.cond_stepwise_causal(memo, context, gen.calculate_pip(memo, p),
                3, .2, z, ld, 100, .05, p0, s, p, allow_dup)
            for j in range(len(start)):
                selected = start[:]
                selected.pop(j)
                gen.add(memo, context, z, ld, 100, .05, p0, s, selected, range(p))
            sets, probabilities, events = gen.find_credible_set(ld, memo, start, p, .95, .05,
                allow_dup, {"purity": 0.}, return_diagnostics=True)
            duration = time.perf_counter() - began
            passes = weights.call_count
        return dict(context=context, memo=memo, pips=gen.calculate_pip(memo, p),
                    null=memo[()], start=start, sets=sets, probabilities=probabilities,
                    events=events, decisions=decisions, seconds=duration, normalizations=passes)

    def assert_equivalent(self, eager, batched):
        self.assertEqual(eager["context"], batched["context"])
        self.assertEqual(list(eager["memo"]), list(batched["memo"]))
        self.assertEqual(eager["memo"], batched["memo"])
        self.assertEqual(eager["null"], batched["null"])
        self.assertEqual(eager["start"], batched["start"])
        self.assertEqual(eager["sets"], batched["sets"])
        self.assertEqual(eager["probabilities"], batched["probabilities"])
        self.assertEqual(eager["events"], batched["events"])
        self.np.testing.assert_array_equal(eager["pips"], batched["pips"])
        self.assertEqual(len(eager["decisions"]), len(batched["decisions"]))
        for a, b in zip(eager["decisions"], batched["decisions"]):
            self.assertEqual(a[:4], b[:4])
            self.np.testing.assert_array_equal(a[4], b[4])
            self.assertEqual(a[5:], b[5:])

    def test_seeded_small_and_hundreds_of_variant_decisions_sets_and_pips_match(self):
        for p, seed, allow_dup in [(12, 31, False), (12, 7, True), (300, 19, False)]:
            with self.subTest(p=p, seed=seed, allow_dup=allow_dup):
                eager = self.expansion(p, seed, True, allow_dup)
                batched = self.expansion(p, seed, False, allow_dup)
                self.assert_equivalent(eager, batched)
                self.assertGreaterEqual(len(batched["start"]), 2,
                    "Equivalence must exercise genuine conditional signal selection")
                self.assertLess(batched["normalizations"], eager["normalizations"])
                record_measurement(dict(scope="real_torch_complete_credible_set_expansion",
                    p=p, seed=seed, allow_duplicates=allow_dup, support_size=len(batched["memo"]),
                    eager_normalizations=eager["normalizations"], batch_normalizations=batched["normalizations"],
                    eager_seconds=eager["seconds"], batch_seconds=batched["seconds"],
                    exact_support_decisions_sets_pips_null=True))

    def test_all_add_delete_change_rounds_match_eager_reference(self):
        gen, pm, torch = self.gen, self.native_pm, self.torch
        p = 8
        z, ld = torch.linspace(-2, 3, p), torch.eye(p)
        s, prior = z.unsqueeze(1), torch.full((p,), .2)
        memo = {(): 0.}
        context = pm.freeze_memo(memo, prior)
        for name in ("add", "delete", "change"):
            eager, ce = copy.deepcopy((memo, context))
            batch, cb = copy.deepcopy((memo, context))
            candidates = []
            real_abf = gen.abf
            def capture(*args, **kwargs):
                candidates.append(tuple(args[7]))
                return real_abf(*args, **kwargs)
            with mock.patch.object(gen, "abf", side_effect=capture):
                getattr(gen, name)(batch, cb, z, ld, 100, .05, prior, s, [0, 2], range(p))
            for selected in candidates:
                real_abf(z, ld, eager, ce, 100, .05, prior, list(selected), s)
            self.assertEqual(cb, ce)
            self.assertEqual(batch, eager)
            self.assertEqual(list(batch), list(eager))

    def test_conditional_membership_and_mask_match_independent_probability_table(self):
        gen, np = self.gen, self.np
        memo = {(): .1, (0,): .2, (1,): .05, (2,): .15,
                (0, 1): .15, (0, 2): .25, (1, 2): .1}
        np.testing.assert_allclose(gen.calculate_pip(memo, 3), [.6, .3, .5], atol=1e-15)
        self.assertEqual(memo[()], .1)
        index, probability, values = gen.find_cond_prob_constrained_causal([0], 2, memo, 3)
        # P(G) = .2+.15+.25=.6; joint masses with 1/2 are .15/.25.
        np.testing.assert_allclose(values, [0., .15/.6, .25/.6], atol=1e-15)
        self.assertEqual(index, 2)
        self.assertAlmostEqual(probability, .25/.6)
        index, probability, masked = gen.find_cond_prob_constrained_causal_no_dup(
            [0], 2, memo, 3, [[2]])
        np.testing.assert_allclose(masked, [0., .15/.6, 0.], atol=1e-15)
        self.assertEqual(index, 1)
        self.assertAlmostEqual(probability, .25)
        # Exclusion masks a member but does not create a prior over remaining choices.
        self.assertLess(float(masked.sum()), 1.)
        self.assertAlmostEqual(gen.conditional_union_probability(memo, [0], [1, 2], 3), 2/3)

    def test_impossible_conditioning_and_malformed_event_inputs_are_rejected(self):
        gen = self.gen
        for method in [lambda: gen.find_cond_prob_constrained_causal([0], 1, {(): 1.}, 2),
                       lambda: gen.find_cond_prob_constrained_causal_no_dup([0], 1, {(): 1.}, 2, []),
                       lambda: gen.conditional_union_probability({(): 1.}, [0], [1], 2)]:
            with self.assertRaisesRegex(ValueError, "positive finite posterior mass"):
                method()
        for given, selected, memo in [([0, 0], [1], {(0,): 1.}),
                ([True], [1], {(0,): 1.}), ([], [2], {(0,): 1.}),
                ([], [0], {(0,): float("nan")}), ([], [0], {(0,): -1.})]:
            with self.assertRaises(ValueError):
                gen.conditional_union_probability(memo, given, selected, 2)

    def test_overlapping_models_use_union_coverage_and_label_legacy_sum(self):
        gen, torch = self.gen, self.torch
        memo = {(0, 1): .4, (0,): .3, (1,): .2, (2,): .1}
        ld = torch.full((3, 3), .2)
        ld.fill_diagonal_(1.)
        sets, probabilities, event = gen.find_credible_set(
            ld, memo, [0], 3, .95, 0., True, {"purity": 0.}, return_diagnostics=True)
        self.assertEqual(sets, [[0, 1, 2]])
        self.assertEqual(probabilities, [[.7, .6, .1]])
        self.assertEqual(event[0]["condition_indices_zero_based"], [])
        self.assertEqual(event[0]["conditional_union_probability"], 1.)
        self.assertTrue(event[0]["target_met"])
        self.assertAlmostEqual(event[0]["legacy_membership_sum"], 1.4)
        sets, _, diagnostic = gen.find_credible_set(ld, memo, [0], 3, .95, 0., True,
            {"purity": 0., "credible_set_rule": "legacy_membership_sum"}, return_diagnostics=True)
        self.assertEqual(sets, [[0, 1]])
        self.assertAlmostEqual(diagnostic[0]["legacy_membership_sum"], 1.3)
        self.assertAlmostEqual(diagnostic[0]["conditional_union_probability"], .9)
        self.assertFalse(diagnostic[0]["target_met"])

    def test_conditional_union_denominator_constraints_and_deduplication_are_honest(self):
        gen, torch = self.gen, self.torch
        memo = {(): .2, (0, 2): .3, (1, 2): .2, (0, 1, 2): .3}
        # Given 2, union{0,1}=1; member marginals are .75/.625, not choices summing1.
        self.assertAlmostEqual(gen.conditional_union_probability(memo, [2], [0], 3), .75)
        self.assertEqual(gen.conditional_union_probability(memo, [2], [0, 1], 3), 1.)
        ld = torch.full((3, 3), .2)
        ld.fill_diagonal_(1.)
        sets, _, diagnostic = gen.find_credible_set(ld, memo, [0, 2], 3, .95, 0., False,
            {"purity": .5}, return_diagnostics=True)
        self.assertEqual(sets[0], [0])
        self.assertEqual(diagnostic[0]["condition_indices_zero_based"], [2])
        self.assertAlmostEqual(diagnostic[0]["conditional_union_probability"], .75)
        self.assertFalse(diagnostic[0]["target_met"])
        out = gen.unique_sets([[0, 1], [1, 0]], [[.6, .4], [.4, .6]],
                              [{"condition": [2]}, {"condition": [3]}])
        self.assertEqual(out, ([[0, 1]], [[.6, .4]], [{"condition": [2]}]))

    def test_main_sidecar_matches_final_frozen_posterior_and_native_variant_order(self):
        gen, pm, np, torch = self.gen, self.native_pm, self.np, self.torch
        z = torch.tensor([4., 3., 2.])
        ld = torch.full((3, 3), .2)
        ld.fill_diagonal_(1.)
        s, prior = torch.inverse(ld) @ z.unsqueeze(1), torch.tensor([.05, .08, .1])
        memo = {(): 0.}
        context = pm.freeze_memo(memo, prior)
        for k in range(1, 4):
            for selected in itertools.combinations(range(3), k):
                gen.abf(z, ld, memo, context, 100, .05, prior, list(selected), s)
        names = ["variant_1", "variant_2", "variant_3"]
        with tempfile.TemporaryDirectory(prefix="fmb-cs-sidecar-") as directory:
            target = Path(directory)
            (target / "data.z").write_text("\n".join(f"{n} {float(x)}" for n, x in zip(names, z)) + "\n")
            np.savetxt(target / "data.ld", ld.numpy())
            with (target / "res").open("wb") as output:
                pickle.dump({"memo_format": pm.MEMO_FORMAT, "posterior_context": context,
                             "memo": memo}, output)
            options = dict(z=str(target/"data.z"), LD=str(target/"data.ld"), target=str(target),
                names=names, sigma_sq=.05, n_sub=100, key_thres=.2, n_causal=2,
                allow_duplicates=False, coverage_ths=.95, selection_prob=0., purity=0.)
            gen.main(options)
            diagnostic = json.loads((target/"credible_set_diagnostics.json").read_text())
            marker = json.loads((ROOT/"BEATRICE_annot_sparse/scripts/credible_set_contract.json").read_text())
            self.assertEqual(diagnostic["format"], marker["format"])
            self.assertEqual(diagnostic["index_base"], marker["index_base"])
            self.assertEqual(diagnostic["stopping_rule"], marker["default_stopping_rule"])
            with (target/"res").open("rb") as source:
                final = pickle.load(source)
            self.assertEqual(diagnostic["variant_names"], names)
            self.assertEqual(diagnostic["index_base"], 0)
            self.assertEqual(diagnostic["posterior_memo_format"], pm.MEMO_FORMAT)
            self.assertEqual(diagnostic["posterior_support_size"], len(final["memo"]))
            self.assertEqual(diagnostic["null_probability"], final["memo"][()])
            self.assertEqual(final["posterior_context"]["final_prior"], context["final_prior"])
            exported_sets = [list(map(int, line.split())) for line in
                             (target/"credible_set.txt").read_text().splitlines()]
            self.assertEqual(exported_sets, [event["set_indices_zero_based"] for event in diagnostic["sets"]])
            for event in diagnostic["sets"]:
                required, members = set(event["condition_indices_zero_based"]), set(event["set_indices_zero_based"])
                denominator = math.fsum(w for c, w in final["memo"].items() if required.issubset(c))
                numerator = math.fsum(w for c, w in final["memo"].items()
                                     if required.issubset(c) and members.intersection(c))
                self.assertEqual(event["conditional_union_probability"], numerator/denominator)
                self.assertEqual(event["target_met"], numerator/denominator >= .95)
            # Keep a genuine main-generated sidecar for the R wrapper's import fixture.
            evidence = os.getenv("FMB_SIDECAR_FIXTURE_DIR")
            if evidence:
                destination = Path(evidence)
                destination.mkdir(parents=True, exist_ok=True)
                for name in ("credible_set_diagnostics.json", "credible_set.txt", "pip.csv"):
                    (destination/name).write_bytes((target/name).read_bytes())
            bad = dict(options, names=list(reversed(names)))
            with self.assertRaisesRegex(ValueError, "identify every variant uniquely"):
                gen.main(bad)


if __name__ == "__main__":
    unittest.main()
