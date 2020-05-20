import argparse
import pickle

import pandas as pd

COMBINATIONS = {
    'all_tissues': ['Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Adrenal_Gland', 'Artery_Aorta',
                    'Artery_Coronary', 'Artery_Tibial', 'Brain_Amygdala', 'Brain_Anterior_cingulate', 'Brain_Caudate',
                    'Brain_Cerebellar', 'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex', 'Brain_Hippocampus',
                    'Brain_Hypothalamus', 'Brain_Nucleus', 'Brain_Putamen', 'Brain_Spinal_cord',
                    'Brain_Substantia_nigra', 'Breast_Mammary_Tissue', 'Cells_Cultured', 'Cells_EBV', 'Colon_Sigmoid',
                    'Colon_Transverse', 'Esophagus_Gastro', 'Esophagus_Mucosa', 'Esophagus_Muscularis', 'Heart_Atrial',
                    'Heart_L_Vent', 'Kidney_Cortex', 'Liver', 'Lung', 'Minor_Salivary', 'Muscle_Skeletal',
                    'Nerve_Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin_Not_Sun_Epsd', 'Skin_Sun_Epsd',
                    'Small_Intestine', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood'],

    'only_brain': ['Brain_Amygdala', 'Brain_Anterior_cingulate', 'Brain_Caudate', 'Brain_Cerebellar',
                   'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex', 'Brain_Hippocampus',
                   'Brain_Hypothalamus', 'Brain_Nucleus', 'Brain_Putamen', 'Brain_Spinal_cord',
                   'Brain_Substantia_nigra'],

    'non_brain': ['Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Adrenal_Gland', 'Artery_Aorta',
                  'Artery_Coronary', 'Artery_Tibial', 'Breast_Mammary_Tissue', 'Cells_Cultured', 'Cells_EBV',
                  'Colon_Sigmoid', 'Colon_Transverse', 'Esophagus_Gastro', 'Esophagus_Mucosa', 'Esophagus_Muscularis',
                  'Heart_Atrial', 'Heart_L_Vent', 'Kidney_Cortex', 'Liver', 'Lung', 'Minor_Salivary', 'Muscle_Skeletal',
                  'Nerve_Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin_Not_Sun_Epsd', 'Skin_Sun_Epsd',
                  'Small_Intestine', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood'],

    'brain_and_wholeblood': ['Brain_Amygdala', 'Brain_Anterior_cingulate', 'Brain_Caudate', 'Brain_Cerebellar',
                             'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex', 'Brain_Hippocampus',
                             'Brain_Hypothalamus', 'Brain_Nucleus', 'Brain_Putamen', 'Brain_Spinal_cord',
                             'Brain_Substantia_nigra', 'Whole_Blood'],

    'brain_and_intestines': ['Brain_Amygdala', 'Brain_Anterior_cingulate', 'Brain_Caudate', 'Brain_Cerebellar',
                             'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex', 'Brain_Hippocampus',
                             'Brain_Hypothalamus', 'Brain_Nucleus', 'Brain_Putamen', 'Brain_Spinal_cord',
                             'Brain_Substantia_nigra', 'Colon_Sigmoid', 'Colon_Transverse', 'Small_Intestine']
}


# This function return the global multilayer matrix for a multilayer matrix
def global_community_multilayer_matrix(partitions, names):
    # Create an empty pandas dataframe with the names of columns and rows given by the genes names
    globalmultiplexity = pd.DataFrame(0, columns=names, index=names)

    for j in range(0, len(names)):
        for k in range(j + 1, len(names)):
            for partition in partitions:
                if names[j] in partition and names[k] in partition:
                    if partition[names[j]] == partition[names[k]]:
                        globalmultiplexity.at[names[j], names[k]] += 1
    return globalmultiplexity


def dictionary_from_multilayer_matrix(gl2, value):
    # Dictionary where I will save the output of the community detection
    mydict = {}
    # In this list I will save the values that I wanna put in the dictionary
    # in this way I can check if I already have a gene and I don't have double genes there
    list_values = []
    for i in range(0, gl2.shape[0]):
        # Create a dictionary where for each gene
        # I have the list of the other genes that are with him tot number of times
        # This will give the name of the row, i.e. of the gene of interest in a certain row that I want to consider
        row_name = gl2.iloc[i, :].name
        # Look for the vector where a certain condition is created
        test = gl2.iloc[i, :] == value
        single_list = []
        # With this I can check if the value that I am considering in the i=th row of the matrix
        # is already there or not
        if row_name not in list_values:
            for j in range(0, gl2.shape[1]):
                if test[j] == True:
                    single_list.append(gl2.columns.values[j])
                    list_values.append(gl2.columns.values[j])
                if len(single_list) > 0:
                    mydict[row_name] = single_list
    return mydict


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("--combination", help="Which combination to try")
    args = parser.parse_args()
    combination = args.combination
    print("Going with", combination)

    # Create an empty list where I will save all the dictionaries of the various genes
    list_dic = []

    # Going over the tissues in that combination
    for tissue_name in COMBINATIONS[combination]:
        print("Reading", tissue_name)

        # Then I can obtain the correlation matrix relative to a given tissue, that we have already obtained before
        # Read the correlation matrix file
        corr_mat = pd.read_pickle("data/corr_" + tissue_name + ".pkl")
        # The columns of the correlation matrix will tell me the name of the genes that are in each community
        genes = corr_mat.columns.values
        # "m" gives the community number/id for each gene.
        # "m" is basically an array with the same length as genes variable and each element has a number
        # corresponding to the community number id
        m, _ = pickle.load(open("results/louvain_modules_" + tissue_name + ".pkl", "rb"))
        # At this point I can create a dictionary to then pass it to the global multiplexity index code
        dictionary = dict(zip(genes, m))
        # Append the dictionary each time at the total list_dic element
        list_dic.append(dictionary)

    # Create a list that will be needed to save the names to give to the global community multilayer matrix
    names_total = []
    for i in range(0, len(list_dic)):
        name = list(list_dic[i].keys())
        names_total = list(set(names_total + name))

    global_multiplexity_matrix = global_community_multilayer_matrix(list_dic, names_total)

    # Save the global multiplexity matrix
    global_multiplexity_matrix.to_csv("results/global_multiplexity_" + combination + "_all.csv")

    # Create a dictionary where I save the genes that always appear together across the various layers
    dict_mat = dictionary_from_multilayer_matrix(global_multiplexity_matrix, len(COMBINATIONS[combination]))

    # WRITE THE FINAL FILE, WHERE EACH ROW IS A SET OF GENES THAT ALWAYS APPEAR TOGETHER
    # Change the name of the file accordingly
    with open("results/global_multiplexity_" + combination + "_genes.csv", 'w') as file:
        for key, val in dict_mat.items():
            file.write(key + "\t" + str(val) + "\n")
