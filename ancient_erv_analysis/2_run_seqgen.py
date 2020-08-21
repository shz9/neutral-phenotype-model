"""
Author: Shadi Zabad
Date: Nov 2018
"""

from Bio import Phylo
import numpy as np
from itertools import islice
import subprocess
import shutil
import os
import errno



def extract_root_seq(paml_dir, root_name):
    """
    This function extracts the root sequence from PAML
    output files.
    """

    f_name = os.path.join(paml_dir, "rst")

    try:
        start_row = subprocess.check_output("grep -n 'List of extant and reconstructed sequences' " + f_name,
                                        shell=True)
        start_row = int(start_row.split(":")[0].strip())
    except Exception:
        raise Exception("Failed to find reconstructed sequences")

    try:
        end_row = subprocess.check_output("grep -n 'Overall accuracy of the' " + f_name,
                                          shell=True)
        end_row = int(end_row.split(":")[0].strip())
    except Exception:
        raise Exception("Failed to find reconstructed sequences")

    root_seq = ""

    with open(f_name, "rb") as fi:
        for row in islice(fi, start_row, end_row):
            if root_name in row:
                root_seq = row.replace(root_name, "").replace(" ", "")
                break

    return root_seq


def create_seqgen_control_files(tree_file,
                                root_seq,
                                branch_scale_factor=1.):

    # ------------------------------
    # Create root_seq control file:

    root_seq_c = str(len(root_seq) - 1) + "\n"
    root_seq_c += "0"*(len(root_seq) - 1) + "\n" # <- all sites are variable and may experience indels.
    root_seq_c += root_seq

    root_seq_file = os.path.join(output_dir, "root.seq")

    with open(root_seq_file, "wb") as outf:
        outf.writelines(root_seq_c)

    # -------------------------------
    # Create tree control file:

    with open(tree_file, "rb") as tf:
        sim_tree = tf.read()

    # -- --
    # Setting simulation parameters:
    max_indel_length = int(.05 * len(root_seq)) # Max indel size is %5 of size of root sequence
    # -- --

    tree_file_content = "[:" + root_seq_file + "]"
    tree_file_content += ("#b" + str(branch_scale_factor) + "#")
    tree_file_content += ("{" + str(max_indel_length) + ",0}") # Use Chang & Benner (2004)
    tree_file_content += sim_tree.strip()
    tree_file_content += "\n"

    tree_output_file = os.path.join(output_dir, "isg.tree")
    with open(tree_output_file, 'wb') as outf:
        outf.writelines(tree_file_content)

    return tree_output_file


def simulate_dna_evolution():

    seqgen_cmds = ['indel-seq-gen',
                   '-m' + model,
                   '-n' + str(num_replicates),
                   '-f' + base_freq,
                   '-e' + output_seq_file,
                   '-w',
                   '-of']

    #if scaling_factor is not None:
    #    seqgen_cmds.append('-s' + scaling_factor)

    seqgen_cmds.append('< ' + '"' + phy_file + '"')

    print " ".join(seqgen_cmds)
    os.system(" ".join(seqgen_cmds))


paml_main_dir = "./data/1_paml_asr/"
input_trees = "./metadata/segment_trees/paml_inferred/"

# Here we simulate under the HKY model only!
num_replicates = 1
model = "HKY"
base_freq = "0.28,0.22,0.22,0.28"

main_output_dir = "./data/2_simulated_sequences/"

# --- --- ---
# Extract the root sequence and perform the 
# forward simulation of DNA sequences:

paml_subdirs = [sd for sd in os.listdir(paml_main_dir)
                if os.path.isdir(os.path.join(paml_main_dir, sd))]

for psd in paml_subdirs:

    paml_input_dir = os.path.join(paml_main_dir, psd)
    output_dir = os.path.join(main_output_dir, psd)
    tree_file = os.path.join(input_trees, psd + ".nwk")

    # --- --- ---
    # Read and transform the tree file:

    tree_obj = Phylo.read(tree_file, "newick")
    root_name = "node #" + str(len(list(tree_obj.get_terminals())) + 1)

    # --- --- ---
    # Create the output directory (if not exists):
    try:
        os.makedirs(output_dir)
    except OSError as ode:
        if ode.errno != errno.EEXIST:
            raise ode

    root_seq = extract_root_seq(paml_input_dir, root_name)
    phy_file = create_seqgen_control_files(tree_file,
                                           root_seq)

    output_seq_file = os.path.join(output_dir, "alignment")

    simulate_dna_evolution()

    os.rename(output_seq_file + ".ma", output_seq_file + ".fa")

