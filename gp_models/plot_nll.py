import sys
from Bio import Phylo
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, "ancient_erv_analysis")
from compute_q_traits import calculate_gc_content
from gp_models import NeutralModel


def plot_nll(segment):
    erv_data = pd.read_csv(f"../ancient_erv_analysis/data/3_q_traits/2_simulated_sequences/{segment}/traits.csv",
                           index_col=0).T
    erv_tree = Phylo.read(f"../ancient_erv_analysis/metadata/segment_trees/paml_inferred/{segment}.nwk", "newick")

    with open(f"../ancient_erv_analysis/data/2_simulated_sequences/{segment}/alignment.root", "r") as algn_f:
        root_seq = str(algn_f.read()).strip()

    true_z0 = calculate_gc_content(root_seq)

    nm = NeutralModel(erv_data['gc_content'], erv_tree,
                      fixed_params={'u': 1.34, 'Zeq': 0.44, 'sigma_eq': 0.2464 / len(root_seq)})
    nm.fit()

    z0 = np.linspace(0., 1., 100)
    L = 800

    psi = (1. / L) * (0.44 + z0 - 2. * z0 * 0.44)

    nll = []

    for i in range(100):
        nm.z0 = z0[i]
        nm.psi = psi[i]
        nll.append(nm.nll())

    plt.plot(z0, nll)
    plt.axvline(true_z0, c='red', ls='--', label='True Z0')
    plt.xlabel("Z0")
    plt.ylabel("NLL")
    plt.title(f"GC Content - {segment}")
    plt.legend()
    plt.show()


plot_nll('segment_320_722')
