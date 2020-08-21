import requests
import os
import io
import json
from multiprocessing import Pool
import pandas as pd
from Bio import Phylo
import itertools


OUTPUT_DIR = "./erv_algn/"
server = "https://rest.ensembl.org"


def download_alignments(chr_table):

    session = requests.Session()

    for ridx, row in chr_table.iterrows():

        reg_desc = "%s:%d-%d" % (row['chr'], row['start'], row['end'])
        print(reg_desc)
        id_ext = "/alignment/region/homo_sapiens/" + reg_desc

        r = session.get(server + id_ext,
                        headers={ "Content-Type" : "application/json"})

        if not r.ok:
            print(r.content)
        else:
            decoded = r.json()[0]
            tr = Phylo.read(io.StringIO(decoded['tree']), "newick")
            species = list(tr.get_terminals())

            if len(species) > 10:
                print("> 10")
                with open(os.path.join(OUTPUT_DIR, reg_desc + ".json"), 'w', encoding='utf-8') as f:
                    json.dump(decoded, f, ensure_ascii=False, indent=4)
            else:
                print(" < 10")


if __name__ == '__main__':

    ann_table = pd.read_csv("Hsap38.txt", sep='\t')
    ann_table = ann_table.loc[(ann_table['N_cnt'] < 1) & (ann_table['AA_length'] > 100), ]

    exclude = ['1', '2', '3', '4']

    chrs = [ann_table.loc[ann_table['chr'] == c, ]
            for c in list(set(ann_table['chr'])) if c not in exclude]

    t_pool = Pool(20)
    inference_results = t_pool.map(download_alignments, chrs)

    t_pool.close()
    t_pool.join()

