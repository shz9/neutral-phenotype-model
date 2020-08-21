from Bio import Phylo
import os
import glob
import re
import io
import pandas as pd
import numpy as np
import itertools
import json
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from multiprocessing import Pool
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def get_clades_from_phylogeny(phy):

    clades = []

    for cl in phy.find_clades(order="postorder"):
        if cl.is_terminal():
            cl.terminals_below = [cl.name]
        else:
            clades.append({"Count_Correction": 1./(len(cl.clades[0].terminals_below)*len(cl.clades[1].terminals_below)),
                           "Pairs": [x + y
                                     for x in cl.clades[0].terminals_below
                                     for y in cl.clades[1].terminals_below]})
            cl.terminals_below = cl.clades[0].terminals_below + cl.clades[1].terminals_below

    return clades


def calculate_frac_a(seq):
    try:
        return float(seq.count("A"))  / (
                seq.count("G") + seq.count("C") +
                seq.count("A") + seq.count("T"))
    except Exception as e:
        return 0.0


def calculate_gc_content(seq):
    try:
        return float(seq.count("G") + seq.count("C")) / (
                seq.count("G") + seq.count("C") +
                seq.count("A") + seq.count("T"))
    except Exception as e:
        print(seq, str(e))
        return 0.0


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


def calculate_distance(seq1, seq2, method='JC'):

    total_sites = 0
    diff_sites = 0

    for i in range(min(len(seq1), len(seq2))):
        if seq1[i] not in ("-", "*", "N") and seq2[i] not in ("-", "*", "N"):
            total_sites += 1
            if seq1[i] != seq2[i]:
                diff_sites += 1

    if total_sites == 0:
        return None

    p = float(diff_sites) / total_sites

    return -0.75*np.log(1. - (4.*p)/3)


def read_parse_json(filename, gap_threshold=0.5):

    with open(filename, 'rb') as f:
        algn = json.load(f)

    tree = Phylo.read(io.StringIO(algn['tree']), 'newick')
    terminal_species = list(tree.get_terminals())

    basename = os.path.basename(filename).replace(".json", "")
    chr_num = basename.split(":")[0]
    start, end = basename.split(":")[1].split("-")

    res_dict = {
        'fname': basename,
        'chr': chr_num,
        'start': int(start),
        'end': int(end),
    }

    br = []
    terminal_species = []

    for clade in tree.find_clades():
        if clade.branch_length is None:
            continue
        br.append(clade.branch_length)
        if clade.is_terminal():
            species_name = "_".join(clade.name.split("_")[:-3])
            for seg in algn['alignments']:
                if seg['species'] in species_name:
                    species_name = seg['species']
            if species_name in terminal_species or species_name in exclude_species:
                tree.prune(clade)
            else:
                clade.name = species_name
                terminal_species.append(species_name)

    # -------------- Tree Statistics --------------
    res_dict['num_terminals'] = len(terminal_species)

    res_dict['min_bl'] = np.min(br)
    res_dict['max_bl'] = np.max(br)
    res_dict['mean_bl'] = np.mean(br)

    pairwise_dist = []

    for sp1, sp2 in itertools.combinations(terminal_species, 2):
        pairwise_dist.append(tree.distance(sp1, sp2))

    res_dict['min_pd'] = np.min(pairwise_dist)
    res_dict['max_pd'] = np.max(pairwise_dist)
    res_dict['mean_pd'] = np.mean(pairwise_dist)

    # -------------- Alignment Statistics --------------

    eff_sp_count = 0
    species_features = {}

    for seg in algn['alignments']:
        seg['species'] = re.sub(r'\[[0-9]*\]', '', seg['species'])
        if seg['species'] in terminal_species:
            if seg['seq'].count("N") == 0 and (float(seg['seq'].count('-')) / len(seg['seq']) <=
                                               gap_threshold):# or len(seg['seq']) -
                                               #seg['seq'].count('-') >= 100):
                eff_sp_count += 1
                species_features[seg['species']] = calculate_gc_content(seg['seq'])

    res_dict['effective_num_terminals'] = eff_sp_count

    res_ddf = []

    for sp1, sp2 in itertools.combinations(species_features.keys(), 2):
        res_ddf.append({
            'Species1': sp1,
            'Species2': sp2,
            'Time': tree.distance(sp1, sp2),
            'Divergence': (species_features[sp1] - species_features[sp2])**2
        })

    res_ddf = pd.DataFrame(res_ddf)

    # -------------- Time-series data  --------------

    paths = []
    root_val = None

    for ts in species_features.keys():

        cum_dist = 0.0
        parent_seq = None

        for cl in tree.get_path(ts):

            try:
                current_seq = [s['seq'] for s in algn['alignments']
                               if cl.name == s['species']][0]
            except Exception as e:
                continue

            if len(current_seq.replace('-', '').replace('N', '')) < 30:
                continue

            if parent_seq is not None:
                dist = calculate_distance(parent_seq, current_seq)

                if dist is None:
                    break

                cum_dist += dist
            else:
                root_val = round(calculate_gc_content(current_seq), 1)

            paths.append({
                'Time': cum_dist,
                'Value': calculate_gc_content(current_seq),
                'Terminal': ts,
                'filename': basename,
                'Root value': root_val
            })

            parent_seq = current_seq

    paths = pd.DataFrame(paths)

    return res_dict, res_ddf, paths


if __name__ == '__main__':

    input_dir = './erv_algn'
    species_tree = "../figure1/metadata/ensemble_species_tree.nwk"
    exclude_species = ['bos_taurus_hybrid', 'sus_scrofa_usmarc', 'bos_indicus_hybrid',
                       'cricetulus_griseus_chok1gshd_scaffold']

    filenames = list(glob.glob(os.path.join(input_dir, "*.json")))

    #read_parse_json(os.path.join(input_dir, "X:146122743-146123246.json"))

    t_pool = Pool(10)
    results = t_pool.map(read_parse_json, filenames)

    t_pool.close()
    t_pool.join()

    res, dfs, paths = [], [], []

    for r in results:
        res.append(r[0])
        if r[0]['effective_num_terminals'] >= 10 and 0.2 <= r[0]['max_pd'] < 3.0:
            dfs.append(r[1])
            paths.append(r[2])


    print(len(dfs))
    res = pd.DataFrame(res)

    div_dfs = pd.concat(dfs)
    #div_dfs = div_dfs.groupby(['Species1', 'Species2'], as_index=False).mean()

    paths = pd.concat(paths)
    paths.to_csv("gc_content_paths.csv")

    unique_species = list(set(list(div_dfs['Species1']) + list(div_dfs['Species2'])))

    gtree = Phylo.read(species_tree, "newick")
    for cl in gtree.get_terminals():
        if cl.name not in unique_species:
            gtree.prune(cl)

    clades = get_clades_from_phylogeny(gtree)

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    sns.scatterplot(x='Time', y='Divergence', data=div_dfs,
                    marker='^', label='data', ax=ax)

    plt.ylim(0, max(div_dfs['Divergence']) + 0.1*max(div_dfs['Divergence']))
    plt.xlabel("Distance")
    plt.ylabel("Divergence")

    plt.savefig("test.pdf")
    plt.close()

    gdf = div_dfs.groupby(['Species1', 'Species2'], as_index=False).mean()
    clade_dist = []
    clade_div = []

    for clade in clades:
        fdf = gdf.loc[(gdf['Species1'] + gdf['Species2']).isin(clade['Pairs']) |
                      (gdf['Species2'] + gdf['Species1']).isin(clade['Pairs']),]

        clade_dist.append(clade['Count_Correction']*fdf['Time'].sum())
        clade_div.append(clade['Count_Correction']*fdf['Divergence'].sum())

    final_df = pd.DataFrame({"Time": clade_dist,
                             "Divergence": clade_div})

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    sns.scatterplot(x='Time', y='Divergence', data=final_df,
                    marker='^', label='data', ax=ax)

    plt.ylim(0, max(final_df['Divergence']) + 0.1*max(final_df['Divergence']))
    plt.xlabel("Distance")
    plt.ylabel("Divergence")

    plt.savefig("clade_test.pdf")
    plt.close()


