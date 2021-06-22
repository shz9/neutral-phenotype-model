"""
Author: Shadi Zabad
Date: October 2018
"""

import numpy as np
import pandas as pd
import subprocess
import glob
import os
import errno


def find_num_replicates(f_name, info_line):

    output = subprocess.check_output("grep -o '" + info_line + "' " + f_name + " | wc -l", shell=True)
    return int(output.strip())


def get_sim_info(f_name):

    with open(f_name, 'r') as dnaf:

        init_line = dnaf.readline()
        td, sl = init_line.strip().split()

    num_rep = find_num_replicates(f_name, init_line.strip())

    return int(td), num_rep


def transform_sequences(fname):

    dirname = os.path.dirname(fname)
    base_name = os.path.basename(fname).replace("_dna_seq.dat", "")
    num_clades, num_rep = get_sim_info(fname)

    with open(os.path.join(dirname, base_name + "_temp_tree.phy"), "r") as f:
        init_seq = f.readlines()[1].split()[1]

    base_freq = [float(bf) for bf in fname.split("/")[-4].split("_")]

    trait_equilibria = {
        trait: func(base_freq)
        for trait, func in trait_equilibria_func.items()
    }

    output_dir = os.path.join(dirname.replace("simulated_sequences", "q_traits"), base_name)

    try:
        os.makedirs(output_dir)
    except OSError as ode:
        if ode.errno != errno.EEXIST:
            raise ode

    for i, rc in enumerate(range(1, num_rep + 1)):
        df = pd.read_csv(fname,
                         sep=r"\s+",
                         skiprows=rc + (rc - 1)*num_clades,
                         nrows=num_clades,
                         header=None,
                         names=['Taxon', 'seq'])

        for trait, func in traits_to_compute.items():

            df[trait] = df['seq'].apply(func)
            df[trait + "_Z0"] = func(init_seq)
            df[trait + "_Zeq"] = trait_equilibria[trait]

        df.to_csv(os.path.join(output_dir, f"{i}.csv"), index=False)


def calculate_equilibrium_gc_content(base_freq):
    return base_freq[1] + base_freq[2]


def calculate_equilibrium_A_percent(base_freq):
    return base_freq[0]


def calculate_equilibrium_gaussian_es(base_freq):
    freq = np.repeat([base_freq], 100, axis=0).T
    return (freq * gauss_effect_sizes).sum()


def calculate_gc_content(seq):
    return float(seq.count("G") + seq.count("C")) / (seq.count("G")
                                                     +
                                                     seq.count("C")
                                                     +
                                                     seq.count("A")
                                                     +
                                                     seq.count("T"))


def calculate_percent_A(seq):
    return float(seq.count("A")) / (seq.count("G")
                                    +
                                    seq.count("C")
                                    +
                                    seq.count("A")
                                    +
                                    seq.count("T"))


def calculate_gaussian_es_trait(seq):

    letters = 'ACGT'
    idx = [letters.index(l) for i, l in enumerate(seq)]

    one_hot = np.zeros(shape=(4, len(seq)))
    one_hot[idx, np.arange(len(seq))] = 1

    return (one_hot * gauss_effect_sizes).sum()


input_dir = "data/simulated_sequences"

gauss_effect_sizes = np.load("metadata/gaussian_effect_sizes.npy")
traits_to_compute = {
    "gc_content": calculate_gc_content,
    "gaussian_es_trait": calculate_gaussian_es_trait,
    "A_percent": calculate_percent_A
}

trait_equilibria_func = {
    "gc_content": calculate_equilibrium_gc_content,
    "gaussian_es_trait": calculate_equilibrium_gaussian_es,
    "A_percent": calculate_equilibrium_A_percent
}


for sim_file in glob.glob(os.path.join(input_dir, "*/*/*/*.dat")):
    print("Computing quantitative traits from:", sim_file)
    transform_sequences(sim_file)
