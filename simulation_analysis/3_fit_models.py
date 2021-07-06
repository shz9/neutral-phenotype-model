import glob
import os
from Bio import Phylo
from gp_models.gp_models import NeutralModel, OU, BM
import pandas as pd
from multiprocessing import Pool


def process_trait_file(trait_f):

    print(trait_f)

    simulation_id = os.path.basename(trait_f).replace(".csv", "")
    _, _, base_freq, tv_ratio, scale, init_seq, _ = trait_f.split("/")
    u = 1. / (1. - sum([float(bf)**2 for bf in base_freq.split("_")]))

    data = pd.read_csv(trait_f, index_col=0)

    if scale == "20":
        tree = Phylo.read("metadata/pruned_mammalian_tree_transformed_20x.nwk", "newick")
    elif scale == "10":
        tree = Phylo.read("metadata/pruned_mammalian_tree_transformed_10x.nwk", "newick")
    elif scale == "5":
        tree = Phylo.read("metadata/pruned_mammalian_tree_transformed_5x.nwk", "newick")
    else:
        tree = Phylo.read("metadata/pruned_mammalian_tree_transformed.nwk", "newick")

    results = []

    for trait in ['gc_content', 'gaussian_es_trait', 'A_percent']:

        true_z0 = data[trait + '_Z0'][0]
        true_zeq = data[trait + '_Zeq'][0]

        nm = NeutralModel(data[trait], tree)
        res = nm.fit()

        inf_z0 = res['Parameters']['Z0']
        inf_zeq = res['Parameters']['Zeq']

        results.append({'Trait': trait,
                        'Simulation': simulation_id,
                        'BaseFreq': base_freq,
                        'TVRatio': tv_ratio,
                        'ScaleFactor': scale,
                        'InitialSequence': init_seq,
                        'Model': 'NeutralModel',
                        'AICc': res['AIC.c'],
                        'True Z0': true_z0,
                        'Inferred Z0': inf_z0,
                        'True Zeq': true_zeq,
                        'Inferred Zeq': inf_zeq,
                        })

        nm = NeutralModel(data[trait], tree, fixed_params={'u': u})
        res = nm.fit()

        inf_z0 = res['Parameters']['Z0']
        inf_zeq = res['Parameters']['Zeq']

        results.append({'Trait': trait,
                        'Simulation': simulation_id,
                        'BaseFreq': base_freq,
                        'TVRatio': tv_ratio,
                        'ScaleFactor': scale,
                        'InitialSequence': init_seq,
                        'Model': f'NeutralModel (fixed u={u:.2f})',
                        'AICc': res['AIC.c'],
                        'True Z0': true_z0,
                        'Inferred Z0': inf_z0,
                        'True Zeq': true_zeq,
                        'Inferred Zeq': inf_zeq,
                        })

        nm = NeutralModel(data[trait], tree, fixed_params={'u': u, 'Zeq': true_zeq})
        res = nm.fit()

        inf_z0 = res['Parameters']['Z0']
        inf_zeq = res['Parameters']['Zeq']

        results.append({'Trait': trait,
                        'Simulation': simulation_id,
                        'BaseFreq': base_freq,
                        'TVRatio': tv_ratio,
                        'ScaleFactor': scale,
                        'InitialSequence': init_seq,
                        'Model': f'NeutralModel (fixed u={u:.2f}, Zeq={true_zeq:.2f})',
                        'AICc': res['AIC.c'],
                        'True Z0': true_z0,
                        'Inferred Z0': inf_z0,
                        'True Zeq': true_zeq,
                        'Inferred Zeq': inf_zeq,
                        })

        ou = OU(data[trait], tree)
        res = ou.fit()

        inf_z0 = res['Parameters']['Z0']
        inf_zeq = res['Parameters']['Zeq']

        results.append({'Trait': trait,
                        'Simulation': simulation_id,
                        'BaseFreq': base_freq,
                        'TVRatio': tv_ratio,
                        'ScaleFactor': scale,
                        'InitialSequence': init_seq,
                        'Model': f'OU',
                        'AICc': res['AIC.c'],
                        'True Z0': true_z0,
                        'Inferred Z0': inf_z0,
                        'True Zeq': true_zeq,
                        'Inferred Zeq': inf_zeq,
                        })

        ou = OU(data[trait], tree, equilibrium_z0=True)
        res = ou.fit()

        inf_z0 = res['Parameters']['Zeq']
        inf_zeq = res['Parameters']['Zeq']

        results.append({'Trait': trait,
                        'Simulation': simulation_id,
                        'BaseFreq': base_freq,
                        'TVRatio': tv_ratio,
                        'ScaleFactor': scale,
                        'InitialSequence': init_seq,
                        'Model': f'OU (Z0=Zeq)',
                        'AICc': res['AIC.c'],
                        'True Z0': true_z0,
                        'Inferred Z0': inf_z0,
                        'True Zeq': true_zeq,
                        'Inferred Zeq': inf_zeq,
                        })

        bm = BM(data[trait], tree)
        bm.fit()

        inf_z0 = res['Parameters']['Zeq']
        inf_zeq = res['Parameters']['Zeq']

        results.append({'Trait': trait,
                        'Simulation': simulation_id, 
                        'BaseFreq': base_freq,
                        'TVRatio': tv_ratio,
                        'ScaleFactor': scale,
                        'InitialSequence': init_seq,
                        'Model': f'BM',
                        'AICc': res['AIC.c'],
                        'True Z0': true_z0,
                        'Inferred Z0': inf_z0,
                        'True Zeq': true_zeq,
                        'Inferred Zeq': inf_zeq,
                        })

    return pd.DataFrame(results)


if __name__ == '__main__':
    pool = Pool(5)

    files_to_process = glob.glob("data/q_traits/*/*/*/*/*.csv")
    print(f"Processing {len(files_to_process)} quantitative measurements...")

    dfs = pool.map(process_trait_file, files_to_process)

    pool.close()
    pool.join()

    res_df = pd.concat(dfs)
    res_df.to_csv("data/inference_results.csv", index=False)
