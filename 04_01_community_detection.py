"""
This script is quite slow and takes a lot of RAM memory.

It is responsible for running the Louvain community algorithm on each tissue. This script needs the pickle files from
the previous script `03`.

All the Louvain communities will be pickled in "results/louvain_modules_*" files, one for each tissue.
"""
import pickle
from multiprocessing import Pool

import bct
import numpy as np
import pandas as pd

from definitions import TISSUES


# python -u 04_01_community_detection.py | tee outputs/output_04_01.txt

def calculate_communitites(tissue_name):
    print(tissue_name + " calculation started...")

    # Creating a dict for all the results to be printed out
    results_dict = {tissue_name: {'q1s': []}}

    corr_mat = pd.read_pickle("data/corr_" + tissue_name + ".pkl")

    corr_mat = corr_mat.replace([np.inf], np.nan).fillna(0)

    corr_arr = corr_mat.values.copy()
    corr_mat = None  # trying to free memory
    corr_arr[(corr_arr > -0.8) & (corr_arr < 0.8)] = 0

    # About negative_asym: https://www.sciencedirect.com/science/article/pii/S105381191100348X
    # Iterative community finetuning.
    m = None
    q0 = -1
    q1 = 0
    while q1 - q0 > 1e-5:
        q0 = q1
        (m, q1) = bct.community_louvain(corr_arr, ci=m, B='negative_asym')
        results_dict[tissue_name]['q1s'].append(q1)
        # print(tissue_name + " q1 = " + str(q1))

    results_dict[tissue_name]['no_com'] = len(np.unique(m))
    # print("Number of communitites: " + str(len(np.unique(m))))

    uniq_arrs = np.unique(np.unique(m, return_counts=True)[1], return_counts=True)
    results_dict[tissue_name]['uniq_arrs'] = uniq_arrs
    # for i in range(len(uniq_arrs[0])):
    #    print("Communitites of size " + str(uniq_arrs[0][i]) + ": " + str(uniq_arrs[1][i]))

    # print("len(np.unique(m)) = " + str(len(np.unique(m))))
    # print("np.unique(np.unique(m, return_counts=True)[1]) = " + str(np.unique(np.unique(m, return_counts=True)[1])))

    pickle.dump((m, q1), open("results/louvain_modules_" + tissue_name + ".pkl", "wb"))

    return results_dict


if __name__ == '__main__':
    NUMBER_THREADS = 2  # This script really takes a lot of RAM so don't make this too high
    pool = Pool(NUMBER_THREADS)

    results = pool.map(calculate_communitites, TISSUES)

    pool.close()
    pool.join()

    # Changing the results to be a dictionary with tissues as keys
    dict_to_print = {}
    for d in results:
        dict_to_print[list(d.keys())[0]] = list(d.values())[0]

    # Printing a final summary for all the tissues
    print("----- FINAL SUMMARY -----")

    for tissue_n in sorted(dict_to_print.keys()):
        print(tissue_n)

        for q1 in dict_to_print[tissue_n]['q1s']:
            print("q1 = " + str(q1))

        print("Number of communitites: " + str(dict_to_print[tissue_n]['no_com']))

        uniq_arrs = dict_to_print[tissue_n]['uniq_arrs']
        for i in range(len(uniq_arrs[0])):
            print("Communitites of size " + str(uniq_arrs[0][i]) + ": " + str(uniq_arrs[1][i]))
