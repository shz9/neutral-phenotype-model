"""
Author: Shadi Zabad

"""

from gp_models import NeutralModel, OU, BM
from Bio import Phylo
import pandas as pd
import glob
import os.path as osp
import sys
sys.path.insert(0, "../ancient_erv_analysis")
from compute_q_traits import (calculate_gc_content, calculate_longest_orf,
                              calculate_percent_A, calculate_longest_orf_alan)

paml_trees = True
simulated_sequences = False

trait_func = {
    "gc_content": calculate_gc_content,
    "longest_ORF": calculate_longest_orf,
    "A_percent": calculate_percent_A,
    "longest_ORF_alan": calculate_longest_orf_alan
}

sim_seq_dir = "../ancient_erv_analysis/data/2_simulated_sequences"

if simulated_sequences:
    traits = "../ancient_erv_analysis/data/3_q_traits/2_simulated_sequences"
else:
    traits = "../ancient_erv_analysis/data/3_q_traits/0_chopped_alignments"

if paml_trees:
    tree_dir = "../ancient_erv_analysis/metadata/segment_trees/paml_inferred"
else:
    tree_dir = "../ancient_erv_analysis/metadata/segment_trees/species"

res = []

for f in glob.glob(osp.join(traits, "*", "traits.csv")):

    print(f)

    seq_name = osp.basename(osp.dirname(f))

    with open(osp.join(sim_seq_dir, seq_name, "alignment.root"), "r") as algn_f:
        root_seq = str(algn_f.read()).strip()

    data = pd.read_csv(f, index_col=0).T
    tree = Phylo.read(osp.join(tree_dir, seq_name + ".nwk"), "newick")

    for t in data.columns:
        print(t)

        if simulated_sequences:
            true_z0 = trait_func[t](root_seq)
        else:
            true_z0 = None

        try:
            m = NeutralModel(data[t], tree)
            fit = m.fit()
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'NeutralModel',
                        'Log-likelihood': fit['Loglikelihood'], 'Corrected AIC': fit['AIC.c'],
                        'LOCO MSE': m.fit_loco()['MSE'][0],
                        'Z0': fit['Parameters']['Z0'], 'True Z0': true_z0, 'Zeq': fit['Parameters']['Zeq']})
        except Exception:
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'NeutralModel',
                        'Log-likelihood': None, 'Corrected AIC': None, 'LOCO MSE': None,
                        'Z0': None, 'True Z0': None, 'Zeq': None})

        try:
            m = NeutralModel(data[t], tree, fixed_params={'u': 1.34})
            fit = m.fit()
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'NeutralModel (fixed u=1.34)',
                        'Log-likelihood': fit['Loglikelihood'], 'Corrected AIC': fit['AIC.c'],
                        'LOCO MSE': m.fit_loco()['MSE'][0],
                        'Z0': fit['Parameters']['Z0'], 'True Z0': true_z0, 'Zeq': fit['Parameters']['Zeq']})
        except Exception:
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'NeutralModel (fixed u=1.34)',
                        'Log-likelihood': None, 'Corrected AIC': None, 'LOCO MSE': None,
                        'Z0': None, 'True Z0': None, 'Zeq': None})

        if t == 'gc_content':
            try:
                m = NeutralModel(data[t], tree, fixed_params={'u': 1.34, 'Z_eq': 0.44})
                fit = m.fit()
                res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'NeutralModel (fixed u=1.34, Z_eq=0.44)',
                            'Log-likelihood': fit['Loglikelihood'], 'Corrected AIC': fit['AIC.c'],
                            'LOCO MSE': m.fit_loco()['MSE'][0],
                            'Z0': fit['Parameters']['Z0'], 'True Z0': true_z0, 'Zeq': fit['Parameters']['Zeq']})
            except Exception:
                res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'NeutralModel (fixed u=1.34, Z_eq=0.44)',
                            'Log-likelihood': None, 'Corrected AIC': None, 'LOCO MSE': None,
                            'Z0': None, 'True Z0': None, 'Zeq': None})

            try:
                m = NeutralModel(data[t], tree, fixed_params={'u': 1.34, 'Z_eq': 0.44,
                                                              'sigma_eq': 0.2464/len(root_seq)})
                fit = m.fit()
                res.append({'Trait': t, 'Sequence': seq_name,
                            'Model': 'NeutralModel (fixed u=1.34, Z_eq=0.44, sigma_eq=0.2464/L)',
                            'Log-likelihood': fit['Loglikelihood'], 'Corrected AIC': fit['AIC.c'],
                            'LOCO MSE': m.fit_loco()['MSE'][0],
                            'Z0': fit['Parameters']['Z0'], 'True Z0': true_z0, 'Zeq': fit['Parameters']['Zeq']})
            except Exception:
                res.append({'Trait': t, 'Sequence': seq_name,
                            'Model': 'NeutralModel (fixed u=1.34, Z_eq=0.44, sigma_eq=0.2464/L)',
                            'Log-likelihood': None, 'Corrected AIC': None, 'LOCO MSE': None,
                            'Z0': None, 'True Z0': None, 'Zeq': None})

        elif t == 'A_percent':
            try:
                m = NeutralModel(data[t], tree, fixed_params={'u': 1.34, 'Z_eq': 0.28})
                fit = m.fit()
                res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'NeutralModel (fixed u=1.34, Z_eq=0.28)',
                            'Log-likelihood': fit['Loglikelihood'], 'Corrected AIC': fit['AIC.c'],
                            'LOCO MSE': m.fit_loco()['MSE'][0],
                            'Z0': fit['Parameters']['Z0'], 'True Z0': true_z0, 'Zeq': fit['Parameters']['Zeq']})
            except Exception:
                res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'NeutralModel (fixed u=1.34, Z_eq=0.28)',
                            'Log-likelihood': None, 'Corrected AIC': None, 'LOCO MSE': None,
                            'Z0': None, 'True Z0': None, 'Zeq': None})

            try:
                m = NeutralModel(data[t], tree, fixed_params={'u': 1.34, 'Z_eq': 0.28,
                                                              'sigma_eq': 0.2016 / len(root_seq)})
                fit = m.fit()
                res.append({'Trait': t, 'Sequence': seq_name,
                            'Model': 'NeutralModel (fixed u=1.34, Z_eq=0.28, sigma_eq=0.2016/L)',
                            'Log-likelihood': fit['Loglikelihood'], 'Corrected AIC': fit['AIC.c'],
                            'LOCO MSE': m.fit_loco()['MSE'][0],
                            'Z0': fit['Parameters']['Z0'], 'True Z0': true_z0, 'Zeq': fit['Parameters']['Zeq']})
            except Exception:
                res.append({'Trait': t, 'Sequence': seq_name,
                            'Model': 'NeutralModel (fixed u=1.34, Z_eq=0.28, sigma_eq=0.2016/L)',
                            'Log-likelihood': None, 'Corrected AIC': None, 'LOCO MSE': None,
                            'Z0': None, 'True Z0': None, 'Zeq': None})

        try:
            m = OU(data[t], tree)
            fit = m.fit()
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'OU',
                        'Log-likelihood': fit['Loglikelihood'], 'Corrected AIC': fit['AIC.c'],
                        'LOCO MSE': m.fit_loco()['MSE'][0],
                        'Z0': fit['Parameters']['Z0'], 'True Z0': true_z0, 'Zeq': fit['Parameters']['Zeq']})
        except Exception:
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'OU',
                        'Log-likelihood': None, 'Corrected AIC': None, 'LOCO MSE': None,
                        'Z0': None, 'True Z0': None, 'Zeq': None})

        try:
            m = OU(data[t], tree, equilibrium_z0=True)
            fit = m.fit()
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'OU (Z0=Zeq)',
                        'Log-likelihood': fit['Loglikelihood'], 'Corrected AIC': fit['AIC.c'],
                        'LOCO MSE': m.fit_loco()['MSE'][0],
                        'Z0': fit['Parameters']['Zeq'], 'True Z0': true_z0, 'Zeq': fit['Parameters']['Zeq']})
        except Exception:
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'OU (Z0=Zeq)',
                        'Log-likelihood': None, 'Corrected AIC': None, 'LOCO MSE': None,
                        'Z0': None, 'True Z0': None, 'Zeq': None})

        try:
            m = BM(data[t], tree)
            fit = m.fit()
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'BM',
                        'Log-likelihood': fit['Loglikelihood'], 'Corrected AIC': fit['AIC.c'],
                        'LOCO MSE': m.fit_loco()['MSE'][0],
                        'Z0': fit['Parameters']['Z0'], 'True Z0': true_z0, 'Zeq': fit['Parameters']['Z0']})
        except Exception:
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'BM',
                        'Log-likelihood': None, 'Corrected AIC': None, 'LOCO MSE': None,
                        'Z0': None, 'True Z0': None, 'Zeq': None})


if simulated_sequences:
    if paml_trees:
        pd.DataFrame(res).to_csv("model_fit_results/simulated_seq_paml_inf_tree.csv")
    else:
        pd.DataFrame(res).to_csv("model_fit_results/simulated_seq_species_tree.csv")
else:
    if paml_trees:
        pd.DataFrame(res).to_csv("model_fit_results/true_seq_paml_inf_tree.csv")
    else:
        pd.DataFrame(res).to_csv("model_fit_results/true_seq_species_tree.csv")

