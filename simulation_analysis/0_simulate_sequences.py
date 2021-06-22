"""
Author: Shadi Zabad
Date: June 2021
"""

import os
import errno


def create_seqgen_seqtree_file(init_seq):

    try:
        os.makedirs(output_dir)
    except OSError as ode:
        if ode.errno != errno.EEXIST:
            raise ode

    with open(tree_file, "rb") as tf:
        sim_tree = tf.read().decode()

    try:
        with open(temp_phy_file, 'w') as outf:
            outf.writelines("1 " + str(len(init_seq)) + "\n")
            outf.writelines("0 " + init_seq + "\n")
            outf.writelines("1\n")
            outf.writelines(sim_tree + "\n")
    except Exception as fce:
        raise fce


def simulate_dna(n, f, t, sf):

    seqgen_cmds = ['seq-gen',
                   '-mHKY',
                   '-k1',
                   '-n' + str(n),
                   '-f' + f,
                   '-s' + str(sf),
                   '-t' + str(t)]

    seqgen_cmds += ['< ' + '"' + temp_phy_file + '"',
                    '> ' + '"' + output_seq_file + '"']

    print(" ".join(seqgen_cmds))
    os.system(" ".join(seqgen_cmds))

# ========================================
# Read and parse commandline arguments

n_replicates = 100
main_output_dir = "./data/simulated_sequences/"
tree_file = "./metadata/pruned_mammalian_tree_transformed.nwk"

initial_sequences = [
    "TGTCAGCCTCCTTTAAGAATTGCCTCCACGAAAGAAGTGCTACCAGGCGGCTCTTCTCCAGGGCAGCTCACACCAGTGACTGATCAACTAGGTATAAAGG",
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG",
    "TAACAAACAACAAAAAGAAATACAAACAAAGAAAAAAGAAAACAGCATCACCAAAAGCAAAGATCAACTAACAGAGCAGCACACACCCAGAGAATAATAA",
    "TGCCACCCCCCTTTACGCCTCCCCCCCCCCGACGACACACAGCTGCCCCCCCCACGGCCCCGCCCTTCTCCCCCGGCACCCCACACCCCGCGCCCCCCCC",
]

frequencies = [
    "0.28,0.22,0.22,0.28",
    "0.25,0.25,0.25,0.25",
    "0.32,0.18,0.18,0.32"
]

tv_ratio = [
    "1.0", "0.5"
]

scaling_factors = [
    "1", "5", "10", "20"
]

# ========================================

for bf in frequencies:
    for tv in tv_ratio:
        for sf in scaling_factors:

            output_dir = os.path.join(main_output_dir, bf.replace(",", "_"), tv, sf)

            for i, init_seq in enumerate(initial_sequences):

                temp_phy_file = os.path.join(output_dir, f'{i}_temp_tree.phy')
                output_seq_file = os.path.join(output_dir, f'{i}_dna_seq.dat')

                create_seqgen_seqtree_file(init_seq)
                simulate_dna(n_replicates, bf, tv, sf)
