"""
This script converts all the gene names from ENSEMBL to GENE IDs. This information is based on the file
`meta_data/genes_ENSEMBL_to_official_gene.csv`, and the data for each tissue comes from the result of script `01_01`
(in data_filtered/ folder).

The resulting CSVs are saved inside data_filtered/ folder.
"""
import numpy as np
import pandas as pd

from definitions import TISSUES

convert_list = pd.read_csv("meta_data/genes_ENSEMBL_to_official_gene.csv", header=None)
convert_list.columns = ["ensemble", "doubt", "geneid"]
convert_list.set_index("ensemble", inplace=True)
convert_list.drop(columns=['doubt'], inplace=True)

old_names = convert_list.index.values
new_names = convert_list['geneid'].values

for tissue_name in TISSUES:
    print("### - " + tissue_name)
    tmp = pd.read_pickle("data_filtered/" + tissue_name + ".pkl")
    print("Columns: " + str(len(tmp.columns.values)))
    print("Unique columns: " + str(len(np.unique(tmp.columns.values))))

    tmp = tmp.rename(columns=lambda x: x.split('.')[0])  # remove numbers after the point
    print("Columns after removing point: " + str(len(tmp.columns.values)))
    print("Unique columns after removing point: " + str(len(np.unique(tmp.columns.values))))

    tmp = tmp.rename(columns=dict(zip(old_names, new_names)))
    print("Columns after translation: " + str(len(tmp.columns.values)))
    print("Unique columns after translation: " + str(len(np.unique(tmp.columns.values))))

    tmp = tmp.loc[:, ~tmp.columns.duplicated()]
    print("Columns after removing duplicates: " + str(len(tmp.columns.values)))
    print("Unique columns after removing duplicates: " + str(len(np.unique(tmp.columns.values))))

    # Removing columns that were not translated
    to_filter = tmp.filter(like="ENSG", axis=1).columns.values
    tmp.drop(columns=to_filter, inplace=True)
    tmp.to_csv("data_filtered/only_geneids_" + tissue_name + ".csv")
