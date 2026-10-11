import torch
from scripts.posterior_memo import (
    add_configuration, configuration_batch, load_context, posterior_weights,
)
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pickle
import pandas as pd
from scripts.convert_to_gpu import gpu
from scripts.convert_to_gpu_and_tensor import gpu_t
from scripts.convert_to_gpu_scalar import gpu_ts
from scripts.convert_to_cpu import cpu
import os
import json
import math
import numbers

matplotlib.use('Agg')


def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, 4)
        
def load_object(filename):
    with open(filename, 'rb') as input:   # Overwrites any existing file.
        obj = pickle.load(input)
    return obj


def cond_prob(M, K, bp):
    total = 0
    for k in M:
        if len(k)==K:
            total += np.squeeze(M[k])
    pip = np.zeros(bp)

    for k in M:
        if len(k)==K:
            for j in k:
                pip[j] += np.squeeze(M[k])
    pip = pip/total
    
    return pip

def find_cond_prob_constrained_causal_no_dup(given, n_caus, M, bp, cred_set):
    """Masked P(j is causal | given are causal); excluded mass is not renormalized.

    ``n_caus`` is a retained legacy argument. The outer stepwise loop controls
    the maximum selected signal count; this function has no cardinality filter.
    """
    
    cred = []
    for c in cred_set:
        cred+=c
    cred = set(cred)    
    total = 0
    ss = given[:]
    # Fractional posterior mass must accumulate in floating point. Integer
    # arrays silently truncated every model contribution to zero.
    p = np.zeros(bp, dtype=float)
    
    for k in M:
        tag = True
        # for cr in cred_set:
        #     if len(np.intersect1d(cr,k))!=1:
        #         tag = False
        if tag:
            if all([s in k for s in ss]):
                total+=M[k]
                for kk in k:
                    if kk not in given and kk not in cred:                   
                        p[kk] += M[k]
        
   
    if not np.isfinite(total) or total <= 0:
        raise ValueError("Conditional probabilities require positive finite posterior mass")
    p = p/total
    
    p = p.squeeze()
    index_max = np.argmax(p)
    
    return index_max, p[index_max],p 

def find_cond_prob_constrained_causal(given, n_caus, M, bp):
    """P(j is causal | given are causal), with given members masked to zero."""
    
    total = 0
    ss = given[:]
    p = np.zeros(bp, dtype=float)
    
    for k in M:
            if all([s in k for s in ss]):
                total+=M[k]
                for kk in k:
                    if kk not in given:
                        p[kk] += M[k]
        
   
    if not np.isfinite(total) or total <= 0:
        raise ValueError("Conditional probabilities require positive finite posterior mass")
    p = p/total
    
    p = p.squeeze()
    index_max = np.argmax(p)
    
    return index_max, p[index_max],p 


def cond_stepwise_causal(M, mean_memo, pip, prior_causal, threshold, Z, LD,\
                         n_sub, sigma_sq, p0, S, bp, allow_dup):
    """ Create key set.
    """
    
    ind       = np.argsort(pip)[::-1]
        
    start_set = sorted(ind[:1])
    
    if pip[start_set[0]]<threshold:
        return []
    
    for nn in range(prior_causal-1):    
        
        add(M, mean_memo, Z, LD, n_sub, sigma_sq, p0, S, sorted(start_set[:]), range(bp))
        
        if allow_dup:
            new_index_snp, index_prob,_ = find_cond_prob_constrained_causal(start_set, nn+2, M, bp)
            
        else:
            
            new_index_snp, index_prob,_ = find_cond_prob_constrained_causal_no_dup(start_set, nn+2, M, bp, [])
            
        if index_prob>threshold:
            start_set.append(new_index_snp)
            #start_set = sorted(start_set)
        else:
            return start_set#sorted(start_set)
    return start_set        #sorted(start_set)
    

                
def _validated_variant_indices(indices, bp):
    values = tuple(indices)
    if (any(not isinstance(j, numbers.Integral) or isinstance(j, bool) for j in values)
            or len(set(values)) != len(values) or any(j < 0 or j >= bp for j in values)):
        raise ValueError("Conditional event indices must be unique integers in range")
    return tuple(int(j) for j in values)


def _conditioned_models(M, given, bp):
    given = _validated_variant_indices(given, bp)
    required = set(given)
    records = []
    for key, weight in M.items():
        key = _validated_variant_indices(key, bp)
        weight = float(weight)
        if not math.isfinite(weight) or weight < 0:
            raise ValueError("Conditional events require finite nonnegative model probabilities")
        if required.issubset(key):
            records.append((frozenset(key), weight))
    denominator = math.fsum(weight for _, weight in records)
    if not math.isfinite(denominator) or denominator <= 0:
        raise ValueError("Conditional probabilities require positive finite posterior mass")
    return given, records, denominator


def _union_probability(records, denominator, members):
    members = set(members)
    return math.fsum(weight for key, weight in records if members.intersection(key)) / denominator


def conditional_union_probability(M, given, members, bp):
    """P(C intersects members | given subset C) within explored posterior support."""
    members = _validated_variant_indices(members, bp)
    _, records, denominator = _conditioned_models(M, given, bp)
    return _union_probability(records, denominator, members)


def find_credible_set(LD, M, start_set, bp, ths, prob_ths, allow_dup, options,
                      return_diagnostics=False):
    """Greedy conditional union-event sets; membership-sum stopping is an ablation.

    Conditioning on the other selected variants and normalization over explored
    configurations do not imply unconditional or full-model-space coverage.
    """
    cred_prob = []
    cred_set =  []
    diagnostics = []
    start_set = list(_validated_variant_indices(start_set, bp))
    rule = options.get('credible_set_rule', 'conditional_union')
    if rule not in ('conditional_union', 'legacy_membership_sum'):
        raise ValueError("Unknown credible-set stopping rule")
    if (not math.isfinite(float(ths)) or not 0 < ths <= 1 or
            not math.isfinite(float(prob_ths)) or not 0 <= prob_ths <= 1 or
            not math.isfinite(float(options['purity'])) or not 0 <= options['purity'] <= 1):
        raise ValueError("Credible-set targets and constraints must be finite probabilities")
    # cred = []
    
    for i, s1 in enumerate(start_set)  :         

        ss = start_set[:]
        ss.pop(i)
        given, records, denominator = _conditioned_models(M, ss, bp)
       
        
        # local cred sets
        cr_set = [s1]
        p_set = []
        
        # posterior cond prob.
        if allow_dup:
            _,_, p = find_cond_prob_constrained_causal(ss, 0, M, bp)
        else:
           
            _,_, p = find_cond_prob_constrained_causal_no_dup(ss, 0, M, bp, cred_set)
        
        ind_cred = np.argsort(p)[::-1]
        
        p_set.append(np.round(p[s1], decimals=3))
        
        # Membership marginals may overlap. Actual event coverage counts each
        # configuration once, even if it contains several candidate members.
        membership_sum = float(p[s1])
        coverage = _union_probability(records, denominator, cr_set)
        
        for it in ind_cred:
                
            if ((rule == 'conditional_union' and coverage >= ths) or
                    (rule == 'legacy_membership_sum' and membership_sum > ths)):
                break
            else:                
                if it != s1 and p[it]>prob_ths:
                    
                    if np.min(np.abs(cpu(LD[cr_set, it]).data.numpy()))> options['purity']:
                        cr_set.append(it)
                        p_set.append(np.round(p[it], decimals=3))
                        membership_sum += float(p[it])
                        coverage = _union_probability(records, denominator, cr_set)
        
        LLD  = LD[:,cr_set]
        LLD = LLD[cr_set,:]
        LLD =np.abs(cpu(LLD).data.numpy())
       # print('Purity = ', LLD)
        #if np.min(LLD)>0.5:        
        
        cred_set.append(cr_set)
        cred_prob.append(p_set)
        diagnostics.append({
            'set_indices_zero_based': [int(j) for j in cr_set],
            'condition_indices_zero_based': list(given),
            'event': 'at_least_one_set_member_is_causal',
            'posterior_conditioning': 'given_variants_are_causal_within_explored_support',
            'conditional_union_probability': float(coverage),
            'target_probability': float(ths),
            'target_met': bool(coverage >= ths),
            'legacy_membership_sum': float(membership_sum),
            'stopping_rule': rule,
            'selection_probability_threshold': float(prob_ths),
            'purity_threshold': float(options['purity']),
        })
    return (cred_set, cred_prob, diagnostics) if return_diagnostics else (cred_set, cred_prob)


def abf( z, ld, memo, mean_memo, n_sub, sigma_sq, p0, ind, S, batch=None):
    '''
    Aproximate bayes factor calculation.

    Parameters
    ----------
    z : z scores.
    ld : LD Matrix.
    memo : Dynamic programming dict.
    mean_memo : Versioned pure-likelihood cache and frozen final prior.
    n_sub : Number of subjects
    sigma_sq : Variance
    p0 : Prior
    ind : indices of causal configurations
    S : Inverse LD multiplied by z.

    '''
    z = z.unsqueeze(1)
    ind_m  = tuple(sorted(ind[:]))
    cc = gpu(torch.ones(len(z)))
    if len(ind)>0:
        if batch is not None and batch.contains(ind_m):
            return None  # a buffered likelihood exists; no probability query
        if batch is None and ind_m in memo:
            return memo[ind_m]
    
        U =  n_sub*torch.diag(sigma_sq*cc)[:,ind]
        V = ld[ind,:]
            
        inv            = torch.inverse(gpu(torch.eye(len(ind))) + torch.mm(V,U))
            
        sigma_inv      = torch.mm(torch.mm(U,inv),V)
            
        sigma          = gpu(torch.eye(len(ind))) + torch.mm(V,U)
            
        sigma2         = torch.matmul(torch.matmul(z.T, sigma_inv),S)/2
        
        # Cache a pure likelihood and reweight the complete search space under
        # the persisted final annotation prior; no uniform-prior fallback.
        log_bf = -torch.logdet(sigma)/2 + sigma2
        if batch is not None:
            batch.add(ind_m, float(log_bf.detach().cpu().numpy().squeeze()))
            return None
        return add_configuration(memo, mean_memo, ind_m,
                                 float(log_bf.detach().cpu().numpy().squeeze()))
    else:
        return None if batch is not None else memo[()]

def delete(memo, mean_memo, z, ld, n_sub, sigma_sq, p0, S, st_Set, search_Set):
    # No selection reads occur inside a complete candidate round. Exit commits
    # all probabilities before the caller makes its next conditional decision.
    with configuration_batch(memo, mean_memo) as batch:
        abf(z, ld, memo, mean_memo, n_sub, sigma_sq, p0, st_Set, S, batch=batch)
        K = list(st_Set)
        for i in range(len(K)):
            kk = K[:]
            kk.pop(i)
            kk = sorted(kk)
            abf(z, ld, memo, mean_memo, n_sub, sigma_sq, p0, kk, S, batch=batch)

def add(memo, mean_memo, z, ld, n_sub, sigma_sq, p0, S, st_Set, search_set):
    with configuration_batch(memo, mean_memo) as batch:
        abf(z, ld, memo, mean_memo, n_sub, sigma_sq, p0, st_Set, S, batch=batch)
        K = list(st_Set)
        for i in search_set:
            kk = K[:]
            if i not in K:
                kk.append(i)
            kk = sorted(kk)
            
            abf(z, ld, memo, mean_memo, n_sub, sigma_sq, p0, kk, S, batch=batch)
    
def change(memo, mean_memo, z, ld, n_sub, sigma_sq, p0, S, st_Set,search_set):
    with configuration_batch(memo, mean_memo) as batch:
        abf(z, ld, memo, mean_memo, n_sub, sigma_sq, p0, st_Set, S, batch=batch)
        K = list(st_Set)
        
        for j, k_con in enumerate(K):
            
            tem_K = K[:]
            tem_K.pop(j)
            for i in search_set:
                kk = tem_K[:]
                if i not in K and i !=k_con:
                    kk.append(i)
                elif i!=k_con:
                    rem_id = kk.index(i)
                    kk.pop(rem_id)
                kk = sorted(kk)
            
                abf(z, ld, memo, mean_memo, n_sub, sigma_sq, p0, kk, S, batch=batch)
    
def calculate_pip(memo,bp):

    pip = np.zeros(bp)
    tot = 0.0
    for k in memo:
        val = float(np.asarray(memo[k]).squeeze())
        tot += val
        for i in k:
            pip[i] += val

    if tot <= 0:
        return pip
    return np.squeeze(pip/tot)


def regularize_ld(LD):
    LD = (LD + LD.T)/2
    s, w = np.linalg.eig(cpu(LD).data.numpy())
    s = np.real(s)
    s_new = torch.zeros(len(s))
    if min(s)<10**-3:
        s_new = torch.ones(len(s))*(min(s)-10**-3)   
        print("\n Adding a constant {} to regularize LD".format(-min(s)+10**-3))
    LD = LD - gpu(torch.diag(s_new))    
    return LD


def unique_sets(cred_set, cred_prob, diagnostics=None):
    D = []
    P = []
    C = []
    E = []
    for i, xx in enumerate(cred_set):
        ids = np.argsort(xx)[::-1]
        x = np.array(xx)[ids]
        if  list(x) not in D:
            D += [list(x)]
            C += [xx]
            P += [cred_prob[i]]
            if diagnostics is not None:
                E.append(diagnostics[i])
    return (C, P, E) if diagnostics is not None else (C, P)
    
    
    x = D

def main(options):
    print('GENERATING CREDIBLE SETS\n')
    ###############################################################################################################################################
    
    z_data = pd.read_table(options['z'], sep=' ', header=None).to_numpy()
    z_variant_names = [str(name) for name in z_data[:, 0]]
    Z = gpu_t(z_data[:, 1].astype(float))
    LD = gpu_t(pd.read_table(options['LD'], sep=' ', header=None).to_numpy())
    LD = regularize_ld(LD)
    
    S = torch.matmul(torch.inverse(LD),Z.unsqueeze(1))
    result = load_object(os.path.join(options['target'], 'res'))
    mean_memo = load_context(result)
    p0 = gpu_t(np.asarray(mean_memo['final_prior'], dtype=float))
    sigma_sq = gpu_ts(options['sigma_sq'])
    bp = len(Z)
    if len(mean_memo['final_prior']) != bp:
        raise ValueError("Final prior does not match the current variant universe")
    m = posterior_weights(mean_memo)

    ###############################################################################################################################################
    pip = calculate_pip(m, bp)
    threshold_causal = options['key_thres']
    prior_n_causal   = min(bp,options['n_causal'] )
    n_sub = options['n_sub']
    allow_dup = options['allow_duplicates']
    start_set = cond_stepwise_causal(m, mean_memo, pip, prior_n_causal, threshold_causal, Z, LD, n_sub, sigma_sq, p0, S, bp, allow_dup)
                    
    ################################################################################################################################################
    
    for ij, ji in enumerate(start_set):
        ss = start_set[:]
        ss.pop(ij)
        add(m, mean_memo, Z, LD, n_sub, sigma_sq, p0, S,ss, range(bp))

    #########################################################################################################################################################
    ths = options['coverage_ths']
    prob_ths = options['selection_prob']
    cred_set, cred_prob, diagnostics = find_credible_set(
        LD, m, start_set, bp, ths, prob_ths, allow_dup, options, return_diagnostics=True)
    
    pip = calculate_pip(m, bp)
    # Expansion changes the explored posterior support. Persist exactly the
    # posterior used by the exported PIPs, under the same frozen learned prior.
    result['memo'] = m
    result['posterior_context'] = mean_memo
    result['mean_memo'] = mean_memo
    result['pip'] = pip
    save_object(result, os.path.join(options['target'], 'res'))
    
    df = {'variant_index':list(range(bp)),'pip':pip, 'variant_names':options['names']}
    df = pd.DataFrame(df)
    df.to_csv(os.path.join(options['target'],'pip.csv'), index=False)
    
    cred_set, cred_prob, diagnostics = unique_sets(cred_set, cred_prob, diagnostics)
    names = options['names']
    if (len(names) != bp or any(not isinstance(name, str) or not name for name in names)
            or len(set(names)) != bp or list(names) != z_variant_names):
        raise ValueError("Credible-set diagnostic names must identify every variant uniquely")
    diagnostic_output = {
        'format': 'beatrice-conditional-credible-set-diagnostics-v1',
        'index_base': 0,
        'variant_count': bp,
        'variant_names': list(names),
        'posterior_memo_format': mean_memo['format'],
        'posterior_support_size': len(m),
        'null_probability': float(m[()]),
        'null_probability_conditioning': 'explored_support_under_frozen_final_prior',
        'stopping_rule': options.get('credible_set_rule', 'conditional_union'),
        'sets': diagnostics,
    }
    with open(os.path.join(options['target'], 'credible_set_diagnostics.json'), 'w') as output:
        json.dump(diagnostic_output, output, allow_nan=False, indent=2)
        output.write('\n')
    cred_str = []
    for item in cred_set:
        s = ''
        for ii in item:
            s += str(ii) + ' '
        cred_str.append(s[:-1])
    
    cred_p = []
    for item in cred_prob:
        s = ''
        for ii in item:
            s += str(ii) + ' '
        cred_p.append(s)
    
    
    f = open(os.path.join(options['target'],'credible_set.txt'),'w')
    f.write('\n'.join(cred_str))
    f.close()
    
    
    f = open(os.path.join(options['target'],'conditional_credible_variants_probability.txt'),'w')
    f.write('\n'.join(cred_p))
    f.close()
    
    C= [np.zeros(3)]
    for i in range(len(cred_set)):
        C.append(np.random.rand(3))
    c =[C[0] for i in range(bp)]
    for ii, i in enumerate(cred_set):
        for j in i:
            c[j] = C[ii+1]
    
    plt.scatter(range(bp), pip, c=c)
    plt.xlabel('variants')
    plt.ylabel('PIP')
    plt.savefig(os.path.join(options['target'],'credible_set.pdf'))
    plt.close()
                
        
        
        
        
        

        
        
        
        
