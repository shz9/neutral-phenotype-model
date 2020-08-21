"""
Author: Shadi Zabad
Date: October 2018
"""

from Bio import AlignIO, Phylo
import os
import errno


input_tree = "./metadata/pruned_mammalian_tree.nwk"
input_alignment = "./data/_ensemble/mafft_mammalian_alignment.fa"

msa_output_dir = "./data/0_chopped_alignments/"
tree_output_dir = "./metadata/segment_trees/species/"

def get_contiguous_segments(msa,
                            max_perc_missing=.2,
                            allowed_gaps=20,
                            min_seg=30,
                            max_stitch_gaps=50,
                            min_seg_len=200):

    """
    * <max_perc_missing>: The maximum proportion of species missing in any given column.
    * <allowed_gaps>: The number of columns we allow gaps (this is is supposed to
        correct for alignment errors and other similar issues).
    * <min_seg>:  Minimum length of contiguous segments to append
    * <max_stitch_gaps>: This parameter is used to stitch contiguous segments.
        If there are 2 contiguous segments that are within
        <max_stitch_gaps> of each other, then stitch them together
    """

    num_species = len(msa)
    seq_len = len(msa[0])

    # The maximum number of species missing in any given column
    max_missing = int(max_perc_missing * num_species)

    segments = []
    start = None

    for i in range(seq_len):
        if msa[:, i].count('-') > max_missing:
            if allowed_gaps > 0:
                allowed_gaps -= 1
            else:
                if start is not None:
                    if i - start >= min_seg:
                        if len(segments) < 1 or start - max_stitch_gaps >= segments[-1][1]:
                            segments.append([start, i])
                        else:
                            segments[-1][1] = i

                allowed_gaps = 20
                start = None

        elif start is None:
            start = i

    if i - max_stitch_gaps > start:
        segments[-1][1] = i

    return [tuple(seg) for seg in segments if seg[1] - seg[0] >= min_seg_len]


def create_segment_msa(msa, start, end,
                       min_per_species=100):

    n_msa = None

    for idx in range(len(msa)):
        if end - start - msa[idx, start:end+1].seq.count("-") >= min_per_species:
            if n_msa is None:
                n_msa = msa[idx:idx+1]
            else:
                n_msa.append(msa[idx])

    n_msa = n_msa[:, start:end+1]

    # Remove empty columns (if any):
    nn_msa = n_msa[:, 0:0]

    for cl in range(len(n_msa[0])):
        if n_msa[:, cl].count("-") != len(n_msa[:, cl]):
             nn_msa += n_msa[:, cl:cl+1]

    seg_output = os.path.join(msa_output_dir, "segment_" + str(start) + "_" + str(end))

    try:
        os.makedirs(seg_output)
    except OSError as ode:
        if ode.errno != errno.EEXIST:
            raise ode

    AlignIO.write(nn_msa,
                  os.path.join(seg_output, "alignment.fa"),
                  "fasta")

    msa_species = [sp.name for sp in nn_msa]

    tree_obj = Phylo.read(input_tree, "newick")

    for tip in tree_obj.get_terminals():
        if tip.name not in msa_species:
            tree_obj.prune(tip)

    Phylo.write(tree_obj,
                os.path.join(tree_output_dir, "segment_" + str(start) + "_" + str(end) + ".nwk"),
                "newick")


try:
    os.makedirs(tree_output_dir)
except OSError as ode:
    if ode.errno != errno.EEXIST:
        raise ode

msa = AlignIO.read(input_alignment, "fasta")

seg_coord = get_contiguous_segments(msa)

for start, end in seg_coord:
    create_segment_msa(msa, start, end)

