"""
Script responsible to create the co-correlation matrices for each tissue, based on the files
data_filtered/only_geneids_CORRECTED_*.csv that came from the `02_01` script.

Each correlation matrix will be pickled inside data/ folder.
"""
from multiprocessing import Pool

import numpy as np
import pandas as pd

from definitions import TISSUES


def calculate_correlation(tissue_name):
    print(tissue_name)
    table = pd.read_csv("data_filtered/only_geneids_CORRECTED_" + tissue_name + ".csv", index_col=0)

    corr_mat = table.corr(method="pearson")

    # Fisher z-transformation
    corr_mat = 0.5 * np.log((1 + corr_mat) / (1 - corr_mat))

    corr_mat.to_pickle("data/corr_" + tissue_name + ".pkl")


if __name__ == "__main__":
    NUMBER_THREADS = 14
    pool = Pool(NUMBER_THREADS)

    pool.map(calculate_correlation, TISSUES)

    pool.close()
    pool.join()
