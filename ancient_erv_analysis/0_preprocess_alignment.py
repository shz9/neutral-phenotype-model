"""
Author: Shadi Zabad
Date: November 2018
"""

from Bio import AlignIO, Phylo, Seq
import subprocess
import os
import errno


input_alignment = "./data/_ensemble/55_eutherian_mammals_alignment.fa"
input_tree = "./metadata/ensemble_species_tree.nwk"

msa = AlignIO.read(input_alignment, "fasta")
tree_obj = Phylo.read(input_tree, "newick")

# Prune the tree to only include species in the alignment:

msa_species = []

for sp in msa:
    sp.id = sp.name = sp.name.split("/")[0]
    sp.description = ""
    msa_species.append(sp.name)

for clade in tree_obj.get_terminals():
    if clade.name not in msa_species:
        tree_obj.prune(clade)

# If there are multiple individuals from the same species or
# closely related species, keep only one and prune the other:

threshold = 0.001 # At least 1 substitution per 1000 base pairs

for clade in tree_obj.find_clades(order='level'):
    if clade.is_preterminal():
        children = sorted([(child, child.branch_length) for child in clade],
                          key=lambda x: x[1])

        if len(children) > 1 and children[0][1] < threshold:
            tree_obj.prune(children[0][0])


tree_tips = [sp.name for sp in tree_obj.get_terminals()]

# Remove pruned species from alignment:
n_msa = None

for idx in range(len(msa)):
    if msa[idx].name in tree_tips and msa[idx].seq.count(".") < .8*len(msa[0]):
        if n_msa is None:
            n_msa = msa[idx:idx+1]
        else:
            n_msa.append(msa[idx])

msa_species = [sp.name for sp in n_msa]

for sp in tree_obj.get_terminals():
    if sp.name not in msa_species:
        tree_obj.prune(sp)

# For MAFFT, replace dots with dashes:

for sp in n_msa:
    sp.seq = Seq.Seq(str(sp.seq).replace(".", "-"))

# Remove empty columns or columns with non-resolved bases:

nn_msa = n_msa[:, 0:0]

for cl in range(len(n_msa[0])):
    if not("N" in n_msa[:, cl] or n_msa[:, cl].count("-") == len(n_msa[:, cl])):
        nn_msa += n_msa[:, cl:cl+1]

# Save pruned tree and filtered alignment:

Phylo.write(tree_obj, "./metadata/pruned_mammalian_tree.nwk", "newick")
AlignIO.write(nn_msa, "./data/_ensemble/filtered_mammalian_alignment.fa", "fasta")

# Run MAFFT:

mafft_call = [
    "mafft",
    "--genafpair",
    "--maxiterate",
    "1000",
    "--preservecase",
    "./data/_ensemble/filtered_mammalian_alignment.fa",
    ">",
    "./data/_ensemble/mafft_mammalian_alignment.fa"
]

os.system(" ".join(mafft_call))

