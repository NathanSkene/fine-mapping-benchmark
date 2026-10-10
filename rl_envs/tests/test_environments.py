import copy
import io
import json
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np

from finemap_rl import TaskFactory, InferenceEnv, ExperimentEnv, MethodDesignEnv
from finemap_rl.tasks import PrivateTask, public_task, write_release, read_release, canonical
from finemap_rl.scoring import brier_loss, panel_loss, validate_prediction
from finemap_rl.baselines import DEFAULT, predict, configuration_posterior, uncertainty_assay, value_information_assay
from finemap_rl.cli import serve, ContextualBandit


class TaskTests(unittest.TestCase):
    def setUp(self): self.factory = TaskFactory(b"unit-test-assessor-key-0001", p=8)

    def test_private_seed_and_split_namespaces(self):
        a = self.factory.make("train", 4); b = self.factory.make("train", 4)
        self.assertEqual(a.inputs, b.inputs); np.testing.assert_array_equal(a.labels, b.labels)
        ids = [self.factory.make(split, 4).group_id for split in ("train", "validation", "test", "design-dev", "design-final")]
        self.assertEqual(len(set(ids)), 5)
        forbidden = {"seed", "scenario", "labels", "beta", "null", "causal_count", "key"}
        self.assertFalse(forbidden & set(a.inputs))
        other = TaskFactory(b"unit-test-assessor-key-0002", p=8).make("train", 4)
        self.assertNotEqual(a.inputs, other.inputs)

    def test_public_schema_rejects_latent_metadata_and_invalid_ld(self):
        x = self.factory.make("train", 2).inputs
        for update in ({"scenario": "null"}, {"seed": 4}, {"truth": [1]}):
            with self.assertRaises(ValueError): public_task({**x, **update})
        bad = copy.deepcopy(x); bad["LD"][0][1] = 5
        with self.assertRaises(ValueError): public_task(bad)
        bad = copy.deepcopy(x); bad["z"][0] = float("nan")
        with self.assertRaises(ValueError): public_task(bad)

    def test_private_label_copy_is_immutable(self):
        t = self.factory.make("train", 1)
        with self.assertRaises(ValueError): t.labels[0] = 1

    def test_release_round_trip_integrity_and_custody(self):
        with tempfile.TemporaryDirectory() as tmp:
            pub, priv = Path(tmp)/"public", Path(tmp)/"private"
            tasks = [self.factory.make("test", i) for i in range(3)]
            write_release(tasks, pub, priv)
            self.assertEqual(read_release(pub), [t.inputs for t in tasks])
            restored = read_release(pub, priv)
            self.assertEqual(restored[0].scenario, tasks[0].scenario)
            self.assertNotIn("labels", (pub/"inputs.json").read_text())
            with self.assertRaises(ValueError): write_release(tasks, pub, priv)
            (pub/"inputs.json").write_text("[]")
            with self.assertRaises(ValueError): read_release(pub)
            with self.assertRaises(ValueError): write_release(tasks, Path(tmp)/"x", Path(tmp)/"x"/"truth")

    def test_annotation_mechanisms_have_declared_direction(self):
        means = {}
        for mechanism in ("annotation_informative", "annotation_reversal"):
            f = TaskFactory(b"unit-test-assessor-key-0001", p=8, scenarios=(mechanism,))
            values = []
            for i in range(200):
                t = f.make("train", i); anno = np.asarray(t.inputs["annotations"])[:, 0]
                values.append(float(anno[t.labels.astype(bool)].mean()))
            means[mechanism] = np.mean(values)
        self.assertGreater(means["annotation_informative"], .4)
        self.assertLess(means["annotation_reversal"], -.4)


class MathematicalTests(unittest.TestCase):
    def setUp(self): self.task = TaskFactory(b"mathematics-key-for-tests", p=8).make("train", 2)

    def test_permutation_and_consistent_sign_recoding(self):
        x = self.task.inputs; original = predict(x)
        order = np.array([7, 2, 0, 4, 6, 1, 3, 5]); signs = np.array([1,-1,1,-1,1,1,-1,1])
        changed = copy.deepcopy(x)
        changed["variant_id"] = [x["variant_id"][i] for i in order]
        changed["z"] = (np.asarray(x["z"])[order]*signs).tolist()
        changed["LD"] = (np.asarray(x["LD"])[np.ix_(order, order)] * np.outer(signs, signs)).tolist()
        changed["annotations"] = np.asarray(x["annotations"])[order].tolist()
        result = predict(changed)
        result_pip, _ = validate_prediction(result, x)
        np.testing.assert_allclose(result_pip, original["pip"], rtol=1e-10, atol=1e-10)
        self.assertAlmostEqual(result["null_probability"], original["null_probability"], places=10)

    def test_perfect_ld_twins_receive_exchangeable_probabilities(self):
        f = TaskFactory(b"mathematics-key-for-tests", p=8, scenarios=("perfect_ld",))
        x = f.make("train", 1).inputs; ld = np.asarray(x["LD"])
        i,j = np.argwhere(np.triu(ld, 1) > 1-1e-10)[0]
        result = predict(x)
        self.assertAlmostEqual(result["pip"][i], result["pip"][j], places=10)

    def test_configuration_mass_and_null_event(self):
        pred, indicators, probabilities = configuration_posterior(self.task.inputs)
        self.assertAlmostEqual(probabilities.sum(), 1.)
        self.assertAlmostEqual(pred["null_probability"], probabilities[0])
        np.testing.assert_allclose(pred["pip"], probabilities @ indicators)
        self.assertAlmostEqual(sum(pred["pip"]), float(probabilities @ indicators.sum(axis=1)))

    def test_positive_weighted_brier_is_proper_in_expectation(self):
        # Analytic Bernoulli expected loss, including an independent null event.
        probabilities = np.array([.15, .6, .3]); null_probability=.2
        def expected(q, r):
            return .5*np.mean(probabilities*(q-1)**2+(1-probabilities)*q**2)+.5*(null_probability*(r-1)**2+(1-null_probability)*r**2)
        honest = expected(probabilities, null_probability)
        for q,r in ((np.zeros(3),1.), (np.ones(3),0.), (probabilities+.1,.3)):
            self.assertGreater(expected(q,r), honest)

    def test_assay_likelihood_updates_odds_in_correct_direction(self):
        x = self.task.inputs; baseline = predict(x)
        att = np.zeros(8); att[0]=1; succ=att.copy()
        yes = predict(x, attempts=att, successes=succ)
        no = predict(x, attempts=att, successes=np.zeros(8))
        self.assertGreater(yes["pip"][0], baseline["pip"][0])
        self.assertLess(no["pip"][0], baseline["pip"][0])

    def test_duplicate_missing_nan_and_cross_task_predictions_fail(self):
        pred = predict(self.task.inputs)
        mutations = [dict(pred, task_id="other"), dict(pred, pip=[float("nan")]*8),
                     dict(pred, variant_id=[pred["variant_id"][0]]*8), dict(pred, pip=[.9]*8,null_probability=.5)]
        for bad in mutations:
            with self.assertRaises(ValueError): brier_loss(bad, self.task)

    def test_failed_tasks_remain_in_denominator(self):
        tasks = [self.task]*3
        report = panel_loss(lambda x: (_ for _ in ()).throw(RuntimeError("failure")), tasks)
        self.assertEqual(report, {"loss": 1., "failures": 3, "attempted": 3})

    def test_missing_generating_target_does_not_force_observed_pip_sum(self):
        t = PrivateTask(self.task.inputs, np.zeros(8), False, "missing_causal", "fixture")
        prediction = {"task_id":t.inputs["task_id"],"variant_id":t.inputs["variant_id"],"pip":[0.]*8,"null_probability":0.}
        self.assertEqual(brier_loss(prediction,t), 0.)


class InteractionTests(unittest.TestCase):
    def setUp(self): self.factory = TaskFactory(b"interaction-key-for-tests", p=8)

    def test_one_shot_no_truth_in_info_and_no_second_score(self):
        env = InferenceEnv(self.factory); obs, info = env.reset(seed=7)
        expected = 1-brier_loss(predict(obs["inputs"]), env._task)
        nxt, reward, done, truncated, info = env.step({"type":"submit","prediction":predict(obs["inputs"])})
        self.assertEqual(reward, expected); self.assertTrue(done); self.assertFalse(truncated)
        self.assertEqual(info,{"status":"submitted"}); self.assertEqual(obs,nxt)
        with self.assertRaises(RuntimeError): env.step({"type":"submit","prediction":predict(obs["inputs"])})

    def test_malformed_action_is_terminal_and_label_independent(self):
        env=InferenceEnv(self.factory); env.reset(seed=1)
        _,reward,done,_,info=env.step({"type":"get_truth"})
        self.assertEqual(reward,-1); self.assertTrue(done); self.assertEqual(info,{"status":"invalid_action"})

    def test_assay_budget_repetition_and_determinism(self):
        env=ExperimentEnv(self.factory,budget=2,max_repeats=1); obs,_=env.reset(seed=9)
        action={"type":"assay","variant_id":obs["inputs"]["variant_id"][0]}
        nxt,reward,done,_,info=env.step(action)
        self.assertEqual(reward,-.003); self.assertFalse(done); self.assertEqual(nxt["budget_remaining"],1)
        self.assertFalse(nxt["legal_assays"][0]); self.assertEqual(sum(nxt["assay_attempts"]),1)
        env.reset(seed=9); repeated=env.step(action)[0]
        self.assertEqual(nxt,repeated)
        _,reward,done,_,_=env.step(action)
        self.assertEqual(reward,-1); self.assertTrue(done)

    def test_active_baseline_terminates_after_budget_plus_submission(self):
        env=ExperimentEnv(self.factory,budget=3); obs,_=env.reset(seed=3); rewards=[]
        for _ in range(4):
            obs,r,done,_,_=env.step(uncertainty_assay(obs)); rewards.append(r)
            if done: break
        self.assertTrue(done); self.assertEqual(len(rewards),4)
        self.assertEqual(rewards[:3],[-.003]*3)

    def test_value_of_information_stops_when_declared_cost_exceeds_gain(self):
        env=ExperimentEnv(self.factory,budget=3,cost=.1); obs,_=env.reset(seed=3)
        action=value_information_assay(obs)
        self.assertEqual(action["type"],"submit")
        # A degenerate declared null prior has zero possible information gain.
        prior={**DEFAULT,"count_prior":[1.,0.,0.]}
        env=ExperimentEnv(self.factory,budget=3,cost=0); obs,_=env.reset(seed=3)
        self.assertEqual(value_information_assay(obs,prior)["type"],"submit")

    def test_design_uses_disjoint_holdout_and_bounded_queries(self):
        env=MethodDesignEnv(self.factory,panel_size=3,query_budget=1); obs,_=env.reset(seed=2)
        self.assertFalse({t.group_id for t in env._dev}&{t.group_id for t in env._final})
        obs,r,done,_,_=env.step({"type":"evaluate","candidate":DEFAULT})
        self.assertEqual(r,-.002); self.assertFalse(done); self.assertEqual(obs["query_budget_remaining"],0)
        self.assertNotIn("labels",json.dumps(obs)); self.assertNotIn("scenario",json.dumps(obs))
        _,reward,done,_,info=env.step({"type":"submit","candidate":DEFAULT})
        self.assertTrue(done); self.assertEqual(info["attempted"],3); self.assertGreaterEqual(reward,0)
        with self.assertRaises(RuntimeError): env.step({"type":"submit","candidate":DEFAULT})

    def test_final_protocol_rejects_chosen_seeds_caps_and_hides_attributes(self):
        env=InferenceEnv(self.factory,"test")
        requests='\n'.join(json.dumps(x) for x in ({"op":"reset","seed":1},{"op":"reset"},{"op":"state"},{"op":"get_truth"},{"op":"reset"}))+"\n"
        output=io.StringIO(); serve(env,io.StringIO(requests),output,split="test",max_episodes=1)
        rows=[json.loads(line) for line in output.getvalue().splitlines()]
        self.assertEqual(rows[0],{"error":"invalid_request"}); self.assertIn("observation",rows[1])
        self.assertEqual(set(rows[2]["state"]),{"active","task_id"})
        self.assertEqual(rows[3],rows[0]); self.assertEqual(rows[4],rows[0])
        self.assertNotIn("labels",output.getvalue()); self.assertNotIn("scenario",output.getvalue())

    def test_nonfinite_json_and_oversized_request_fail_closed(self):
        env=InferenceEnv(self.factory); output=io.StringIO()
        serve(env,io.StringIO('{"op":"reset","seed":NaN}\n'),output)
        self.assertEqual(json.loads(output.getvalue()),{"error":"invalid_request"})
        output=io.StringIO(); serve(env,io.StringIO('x'*1_048_580+'\n'),output)
        self.assertEqual(len(output.getvalue().splitlines()),1)

    def test_final_task_cannot_be_skipped_by_resetting_while_active(self):
        env=InferenceEnv(self.factory,"test"); output=io.StringIO()
        serve(env,io.StringIO('{"op":"reset"}\n{"op":"reset"}\n{"op":"state"}\n'),output,split="test",max_episodes=5)
        rows=[json.loads(x) for x in output.getvalue().splitlines()]
        self.assertEqual(rows[1],{"error":"invalid_request"})
        self.assertEqual(rows[0]["observation"]["inputs"]["task_id"],rows[2]["state"]["task_id"])

    def test_bandit_parameters_are_updated_from_actual_reward(self):
        learner=ContextualBandit([DEFAULT]); task=self.factory.make("train",3)
        self.assertEqual(learner.choose(task.inputs),0)
        learner.update(task.inputs,0,.9)
        self.assertEqual(learner.counts.tolist(),[1]); self.assertGreater(np.linalg.norm(learner.b),0)


class GymTests(unittest.TestCase):
    def test_official_checker_all_numeric_adapters(self):
        try:
            from gymnasium.utils.env_checker import check_env
            from finemap_rl.gym_adapter import FineMappingGym, MethodDesignGym
        except ImportError: self.skipTest("Install the gym extra to verify compatibility")
        for env in (FineMappingGym(p=8), FineMappingGym(p=8,experiment=True),MethodDesignGym(p=8,panel_size=3,query_budget=2)):
            check_env(env,skip_render_check=True)
            env.close()


class OpenEnvTests(unittest.TestCase):
    def test_typed_adapter_and_persistent_websocket(self):
        try:
            from fastapi.testclient import TestClient
            from finemap_rl.openenv_adapter import create_finemap_app
        except ImportError: self.skipTest("Install the openenv extra for hosted API compatibility")
        with TestClient(create_finemap_app(p=8,split="test",max_episodes=1)) as client:
            self.assertEqual(client.get("/health").status_code,200)
            self.assertEqual(client.get("/schema").status_code,200)
            with client.websocket_connect("/ws") as ws:
                ws.send_json({"type":"reset","data":{"seed":2}})
                self.assertEqual(ws.receive_json()["type"],"error")
                ws.send_json({"type":"reset","data":{}})
                data=ws.receive_json()["data"]
                x=data["observation"]["payload"]["inputs"]
                ws.send_json({"type":"state"})
                state=ws.receive_json()["data"]
                self.assertEqual(set(state),{"episode_id","step_count","active"})
                self.assertTrue(state["active"])
                ws.send_json({"type":"step","data":{"command":{"type":"submit","prediction":predict(x)}}})
                data=ws.receive_json()["data"]
                self.assertTrue(data["done"]); self.assertGreaterEqual(data["reward"],0)
                self.assertNotIn("labels",json.dumps(data)); self.assertNotIn("scenario",json.dumps(data))
                ws.send_json({"type":"step","data":{"command":{"type":"submit","prediction":predict(x)}}})
                self.assertEqual(ws.receive_json()["type"],"error")
                ws.send_json({"type":"reset","data":{}})
                self.assertEqual(ws.receive_json()["type"],"error")

    def test_sessions_use_independent_private_keys(self):
        try:
            from fastapi.testclient import TestClient
            from finemap_rl.openenv_adapter import create_finemap_app
        except ImportError: self.skipTest("Install the openenv extra")
        with TestClient(create_finemap_app(p=8)) as client:
            with client.websocket_connect("/ws") as first, client.websocket_connect("/ws") as second:
                first.send_json({"type":"reset","data":{"seed":3}})
                second.send_json({"type":"reset","data":{"seed":3}})
                a=first.receive_json()["data"]["observation"]["payload"]["inputs"]
                b=second.receive_json()["data"]["observation"]["payload"]["inputs"]
                self.assertNotEqual(a["task_id"],b["task_id"])


@unittest.skipUnless(os.environ.get("FINEMAP_RL_DOCKER")=="1", "Opt-in actual Docker integration/red-team tests")
class DockerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from finemap_rl.code_runner import DockerPredictorRunner
        cls.runner=DockerPredictorRunner(timeout=4)
        cls.task=TaskFactory(b"docker-test-private-key-00",p=8).make("train",1)

    def test_real_source_submission_and_heldout_scoring(self):
        from finemap_rl.code_runner import STARTER_CODE,CodeImprovementEnv
        outputs=self.runner.run(STARTER_CODE,[self.task.inputs])
        validate_prediction(outputs[0],self.task.inputs)
        env=CodeImprovementEnv(self.runner,factory=TaskFactory(b"docker-test-private-key-00",p=8),panel_size=3,query_budget=1)
        env.reset(seed=1); obs,r,done,_,_=env.step({"type":"evaluate","source":STARTER_CODE})
        self.assertFalse(done); self.assertEqual(r,-.002); self.assertEqual(obs["validation_history"][0]["metrics"]["failures"],0)
        _,r,done,_,info=env.step({"type":"submit","source":STARTER_CODE})
        self.assertTrue(done); self.assertEqual(info["failures"],0); self.assertGreaterEqual(r,0)

    def test_no_assessor_files_or_network_or_writable_source(self):
        from finemap_rl.code_runner import STARTER_CODE
        probes='''import os, socket
assert os.getuid() == 65534
assert not os.path.exists('/private/truth.json')
assert not os.path.exists('/var/run/docker.sock')
assert 'ASSESSOR_KEY' not in os.environ
try:
    open('/agent/predictor.py','w').write('bad')
    raise AssertionError('source writable')
except OSError: pass
s=socket.socket(); s.settimeout(.5)
try:
    s.connect(('1.1.1.1',80))
    raise AssertionError('network available')
except OSError: pass
'''
        outputs=self.runner.run(probes+STARTER_CODE,[self.task.inputs])
        validate_prediction(outputs[0],self.task.inputs)

    def test_timeout_output_flood_invalid_json_and_hidden_file_attempts(self):
        for source in ('while True: pass', 'import os\nwhile True: os.write(1,b"x"*4096)',
                       'print("not json")', "open('/private/truth.json').read()"):
            with self.assertRaises((RuntimeError,ValueError)):
                self.runner.run(source,[self.task.inputs])

    def test_failed_code_charged_for_every_attempted_locus(self):
        from finemap_rl.code_runner import CodeImprovementEnv
        env=CodeImprovementEnv(self.runner,factory=TaskFactory(b"docker-test-private-key-00",p=8),panel_size=3)
        env.reset(seed=4); obs,_,_,_,_=env.step({"type":"evaluate","source":"raise RuntimeError('bad')"})
        report=obs["validation_history"][0]["metrics"]
        self.assertEqual(report,{"loss":1.,"failures":3,"attempted":3})
        _,reward,done,_,info=env.step({"type":"submit","source":"print('bad')"})
        self.assertEqual(reward,0); self.assertTrue(done); self.assertEqual(info["failures"],3)


if __name__ == "__main__": unittest.main()
