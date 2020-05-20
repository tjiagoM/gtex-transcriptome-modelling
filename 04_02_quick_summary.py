import pickle

import numpy as np
import pandas as pd

from definitions import TISSUES

# python -u 04_02_quick_summary.py | tee outputs/output_04_02.txt

# Just to check whether communities are disconnected
for tissue_name in TISSUES:
    print("## " + tissue_name)

    communities, _ = pickle.load(open("results/louvain_modules_" + tissue_name + ".pkl", "rb"))

    corr_mat = pd.read_pickle("data/corr_" + tissue_name + ".pkl")
    corr_mat = corr_mat.replace([np.inf], np.nan).fillna(0)
    corr_arr = corr_mat.values.copy()
    corr_arr[(corr_arr > -0.8) & (corr_arr < 0.8)] = 0

    size_connections = {}

    uniqs = np.unique(communities, return_counts=True)
    for ind, size in enumerate(uniqs[1]):
        if size == 1:
            id_comunity = uniqs[0][ind]

            pos = np.where(communities == id_comunity)

            sum_val = np.sum(corr_arr[pos, :])
            # If it is connected to other communities, it will print
            if sum_val > 0:
                print("Community of size 1 connected! ... " + str(sum_val))

        # Checking whether communities are disconnected among them
        else:
            # Community IDs of a certain size
            com_sizes = [uniqs[0][i] for i, elem in enumerate(uniqs[1]) if elem == size]
            com_others = [uniqs[0][i] for i, elem in enumerate(uniqs[1]) if elem != size and elem != 1]

            pos_sizes = np.where(np.isin(communities, com_sizes))
            pos_others = np.where(np.isin(communities, com_others))

            all_df = pd.DataFrame(corr_arr, index=corr_mat.index, columns=corr_mat.columns)

            filt_indexes = all_df.index.values[pos_sizes]
            filt_columns = all_df.index.values[pos_others]

            # DataFrame where indexes are genes in communities of size `size`
            #   and columns are all the other genes (except of belonging to communities of size 1)
            filt_df = all_df.loc[filt_indexes, filt_columns]

            sumed_df = filt_df.sum(axis=1)
            sum_shape = sumed_df[sumed_df != 0].shape

            size_connections[size] = sum_shape[0]

    for size, shape in sorted(size_connections.items()):
        print("There are " + str(shape) + " connections to outside the communities of size " + str(size))
