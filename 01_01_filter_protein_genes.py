"""
Script responsible to remove genes from the data which do not encode a protein.

This is based on information in `meta_data/uniq_protein_encoding_genes.txt` and saved inside data_filtered/ folder.
"""
import pandas as pd

from definitions import FILES

text_file = open("meta_data/uniq_protein_encoding_genes.txt", "r")
protein_genes = text_file.read().split('\n')
text_file.close()

protein_genes = [x.split('.')[0] for x in protein_genes]

for f in sorted(FILES.items()):
    tmp = pd.read_csv(f[1], sep="\t").drop(columns=["#chr", "start", "end"]).set_index("gene_id").transpose()
    # Getting the genes from this tissue that encode proteins, and filter only those genes
    common = [x for x in tmp.columns.values if x.split('.')[0] in protein_genes]

    filtered_df = tmp.loc[:, common]

    filtered_df.to_pickle("data_filtered/" + f[0] + ".pkl")
