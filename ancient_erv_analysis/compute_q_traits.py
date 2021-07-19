"""
Author: Shadi Zabad
Date: October 2018
"""

import pandas as pd
from Bio import AlignIO, Phylo
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import glob
import shutil
import os
import errno


input_dirs = [
    "./data/0_chopped_alignments/",
    "./data/2_simulated_sequences/"
]

tree_dir = "./metadata/segment_trees/species/"
output_dir_main = "./data/3_q_traits/"


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


def calculate_longest_orf(seq):
    """
    Parts of this function are borrowed from here:
    http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc302
    """

    seq = Seq(str(seq).replace("-", ""), generic_dna)
    max_len = 0

    for strand, nuc in [(1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(seq)-frame) // 3)
            max_len = max(max_len, max([len(pro) for pro in
                                        nuc[frame:frame+length].translate().split("*")]))

    return float(max_len) / len(seq)


def calculate_longest_orf_alan(seq):

    s = str(seq).replace("-", "")

    codons = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S",
              "TCC": "S", "TCA": "S", "TCG": "S", "TAT": "Y", "TAC": "Y",
              "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*",
              "TGG": "W", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
              "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAT": "H",
              "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R",
              "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I",
              "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
              "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S",
              "AGC": "S", "AGA": "R", "AGG": "R", "GTT": "V", "GTC": "V",
              "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A",
              "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
              "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G", "---": "-"
              }

    def revcom(seq):
        com = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'}
        return ''.join(reversed([com[b] for b in list(seq)]))

    maxpep = ""
    for seq in [s, revcom(s)]:
        bases = list(seq)
        for frame in [0, 1, 2]:
            thisp = []
            cDNA = ""
            for i in range(frame, len(bases) - 3, 3):
                trip = ''.join(bases[i:i + 3])

                aa = codons[trip]
                if aa == '*':
                    cDNA += trip
                    # if (len(thisp)>10): print (i-len(cDNA),i,cDNA)
                    if len(thisp) > len(maxpep): maxpep = "".join(thisp)
                    thisp = []
                    cDNA = ""
                elif len(thisp) > 0:
                    thisp.append(aa)
                    cDNA += trip
                elif aa == 'M':
                    thisp.append(aa)
                    cDNA += trip
            if len(thisp) > len(maxpep): maxpep = "".join(thisp)  ## in case no stop codon

    return len(maxpep)


traits_to_compute = [
    ("gc_content", calculate_gc_content),
    ("longest_ORF", calculate_longest_orf),
    ("A_percent", calculate_percent_A),
    ("longest_ORF_alan", calculate_longest_orf_alan)
]


for input_dir in input_dirs:

    parent_dir = os.path.join(output_dir_main,
                              os.path.basename(os.path.normpath(input_dir)))

    for subdir in glob.glob(os.path.join(input_dir, "*/")):

        algn_name = os.path.basename(os.path.normpath(subdir))

        output_dir = os.path.join(parent_dir, algn_name)

        try:
            os.makedirs(output_dir)
        except OSError as ode:
            if ode.errno != errno.EEXIST:
                raise ode

        tree_obj = Phylo.read(os.path.join(tree_dir,
                                           algn_name + ".nwk"), "newick")
        tips = [tip.name for tip in tree_obj.get_terminals()]

        alignment = AlignIO.read(os.path.join(subdir, "alignment.fa"),
                                 "fasta")
        species_traits = {}

        for al in alignment:
            if al.name in tips:
                species_traits[al.name] = [func(al.seq) for trait_name, func in traits_to_compute]

        df = pd.DataFrame(species_traits)
        df.index = [t[0] for t in traits_to_compute]

        df.to_csv(os.path.join(output_dir, "traits.csv"))
