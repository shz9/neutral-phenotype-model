from gp_models import NeutralModel, OU, BM
from Bio import Phylo
import pandas as pd
import glob
import os.path as osp

traits = "../ancient_erv_analysis/data/3_q_traits/2_simulated_sequences"
tree_dir = "../ancient_erv_analysis/metadata/segment_trees/paml_inferred"

res = []

for f in glob.glob(osp.join(traits, "*", "traits.csv")):

    print(f)

    seq_name = osp.basename(osp.dirname(f))

    data = pd.read_csv(f, index_col=0).T
    tree = Phylo.read(osp.join(tree_dir, seq_name + ".nwk"), "newick")

    for t in data.columns:
        print(t)

        try:
            m = NeutralModel(data[t], tree)
            fit = m.fit()
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'NeutralModel',
                        'Corrected AIC': fit['AIC.c'], 'LOCO MSE': m.fit_loco()['MSE'][0],
                        'Z0': fit['Parameters']['Z0'], 'Zeq': fit['Parameters']['Zeq']})
        except Exception:
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'NeutralModel',
                        'Corrected AIC': None, 'LOCO MSE': None, 'Z0': None, 'Zeq': None})

        try:
            m = NeutralModel(data[t], tree, u=1.3)
            fit = m.fit()
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'NeutralModel (fixed u=1.3)',
                        'Corrected AIC': fit['AIC.c'], 'LOCO MSE': m.fit_loco()['MSE'][0],
                        'Z0': fit['Parameters']['Z0'], 'Zeq': fit['Parameters']['Zeq']})
        except Exception:
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'NeutralModel (fixed u=1.3)',
                        'Corrected AIC': None, 'LOCO MSE': None,
                        'Z0': None, 'Zeq': None})

        try:
            m = OU(data[t], tree)
            fit = m.fit()
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'OU',
                        'Corrected AIC': fit['AIC.c'], 'LOCO MSE': m.fit_loco()['MSE'][0],
                        'Z0': fit['Parameters']['Z0'], 'Zeq': fit['Parameters']['Zeq']})
        except Exception:
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'OU',
                        'Corrected AIC': None, 'LOCO MSE': None,
                        'Z0': None, 'Zeq': None})

        try:
            m = OU(data[t], tree, equilibrium_z0=True)
            fit = m.fit()
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'OU (Z0=Zeq)',
                        'Corrected AIC': fit['AIC.c'], 'LOCO MSE': m.fit_loco()['MSE'][0],
                        'Z0': fit['Parameters']['Zeq'], 'Zeq': fit['Parameters']['Zeq']})
        except Exception:
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'OU (Z0=Zeq)',
                        'Corrected AIC': None, 'LOCO MSE': None,
                        'Z0': None, 'Zeq': None})

        try:
            m = BM(data[t], tree)
            fit = m.fit()
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'BM',
                        'Corrected AIC': fit['AIC.c'], 'LOCO MSE': m.fit_loco()['MSE'][0],
                        'Z0': fit['Parameters']['Z0'], 'Zeq': fit['Parameters']['Z0']})
        except Exception:
            res.append({'Trait': t, 'Sequence': seq_name, 'Model': 'BM',
                        'Corrected AIC': None, 'LOCO MSE': None,
                        'Z0': None, 'Zeq': None})


pd.DataFrame(res).to_csv("simulated_seq_paml_inf_tree.csv")

