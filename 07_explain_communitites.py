"""
Script responsible to enrich all the communities (with size greater than 3) using the gseapy package.

Results are saved inside results/EnrichClass/ folder for each community and inside results/enrichment/ for the
significant enrichments.
"""
import argparse
import pickle
from time import sleep

import gseapy as gp
import numpy as np
import pandas as pd

from definitions import TISSUES


def enrich_tissue(tissue_name, last_num):
    communities = pickle.load(open("results/louvain_modules_" + tissue_name + ".pkl", "rb"))
    corr_mat = pd.read_pickle("data/corr_" + tissue_name + ".pkl")

    community_id = 1

    for community in np.unique(communities[0]):
        common = np.array(corr_mat.columns)[communities[0] == community]

        if len(common) <= 3:
            continue
        # print(community_id)
        print("For community", community_id, "(len: " + str(len(common)) + ")...")

        if community_id < last_num:
            community_id += 1
            continue

        enr = gp.enrichr(gene_list=list(common.astype('<U3')),
                         organism='human',
                         description=tissue_name + "_" + str(community_id),
                         gene_sets='Reactome_2016',
                         cutoff=0.05,
                         outdir='results/EnrichClass')
        if enr.results.shape[0] > 0:
            enr.results = enr.results[enr.results['Adjusted P-value'] < 0.05]
            if enr.results.shape[0] > 0:
                enr.results.to_csv(
                    "results/enrichment/" + tissue_name + "_" + str(community_id) + "_" + str(len(common)) + ".csv")
                print("Enriched!")

        community_id += 1
        sleep(50)  # just to go easy on the Enrich API... (constantly getting errors after a while)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--tissue_num", type=int,
                        help='Tissue number in the TISSUES array (from definitions module), on which the code will be executed')
    parser.add_argument("--last_num", type=int,
                        help='Community number from where to start (eg. in case previous run ended unexpectedly')
    args = parser.parse_args()
    print("Going with", TISSUES[args.tissue_num])

    enrich_tissue(TISSUES[args.tissue_num], args.last_num)
