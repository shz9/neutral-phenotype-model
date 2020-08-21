"""
Author: Shadi Zabad
Date: October 2018
"""

from Bio import Phylo
from Bio.Phylo.PAML import baseml
import glob
import json
import os
import errno


sample_ctrl = "./metadata/baseml.ctl"
input_dir = "./data/0_chopped_alignments/"
tree_input_dir = "./metadata/segment_trees/species/"

output_dir = "./data/1_paml_asr/"
tree_output_dir = "./metadata/segment_trees/paml_inferred/"


def run_baseml(alignment_dir):

    algn_name = os.path.basename(os.path.normpath(alignment_dir))
    alignment_file = os.path.join(alignment_dir, "alignment.fa")
    tree_file = os.path.join(tree_input_dir, algn_name + ".nwk")

    outdir = os.path.join(output_dir,
                          algn_name)

    try:
        os.makedirs(outdir)
    except OSError as ode:
        if ode.errno != errno.EEXIST:
            raise ode

    tree_obj = Phylo.read(tree_file, "newick")
    mash_tree_dist = 0.0

    for cl in tree_obj.find_clades():
        mash_tree_dist += cl.branch_length

    bml = baseml.Baseml(alignment=alignment_file, tree=tree_file,
                        working_dir=outdir)

    bml.read_ctl_file(sample_ctrl)
    bml.set_options(model=4)
    bml.out_file = os.path.join(outdir, "result.txt")

    res = bml.run()

    with open(os.path.join(tree_output_dir, algn_name + ".nwk"), "wb") as tf:
        tf.write(res['tree'].replace(" ", ""))

    diff_dict = {
        "Species Distance (Mash)": mash_tree_dist,
        "Region Distance": res['tree length'],
        "Ratio": res['tree length'] / mash_tree_dist
    }

    with open(os.path.join(outdir, "distance_stats.json"), "wb") as dfs:
        json.dump(diff_dict, dfs)


try:
    os.makedirs(tree_output_dir)
except OSError as ode:
    if ode.errno != errno.EEXIST:
        raise ode

for a_dir in glob.glob(input_dir + "*/"):
    run_baseml(a_dir)

