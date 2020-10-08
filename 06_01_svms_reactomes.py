"""
This script runs the SVM prediction analysis from the paper for each REACTOME pathway, which is passed as a flag.

Results are pickled in the svm_results/ folder.
"""
import argparse
import pickle

import pandas as pd
from sklearn import svm
from sklearn.model_selection import cross_validate

from definitions import TISSUES


# python -u 06_01_svms_reactomes.py --react_name NAME | tee outputs/react/output_06_01_NAME.txt


def calculate_react_svm(react):
    reactome_genes = pd.read_table("meta_data/ReactomeData/" + react + ".txt", header=None)

    # Getting only the reactome genes that are actually present in all_df
    common = [x for x in all_df.columns.values if x in reactome_genes.iloc[:, 0].values]
    filtered_df = all_df.loc[:, common].copy()

    filtered_w_tissue = filtered_df.join(all_df['tissue'])

    dic_community = {}

    X = filtered_w_tissue.loc[:, common].values

    for f_n in TISSUES:
        print("-- Predicting " + f_n + " with " + react)
        y = filtered_w_tissue.loc[:, 'tissue'].copy().values
        for j, elem in enumerate(y):
            if elem == "_" + f_n:
                y[j] = 1
            else:
                y[j] = 0

        dic_community[f_n] = {}

        scoring = ['accuracy', 'f1', 'roc_auc']

        clf = svm.SVC(kernel='linear', class_weight="balanced")
        scores = cross_validate(clf, X, list(y), cv=3, scoring=scoring, n_jobs=-1)

        score = scores['test_accuracy']
        dic_community[f_n]['acc'] = score.mean()
        dic_community[f_n]['acc_std'] = score.std()
        print("Accuracy: %.4f (%.4f)" % (score.mean(), score.std()))

        score = scores['test_f1']
        dic_community[f_n]['f1'] = score.mean()
        dic_community[f_n]['f1_std'] = score.std()
        print("F1 score: %.4f (%.4f)" % (score.mean(), score.std()))

        score = scores['test_roc_auc']
        dic_community[f_n]['roc'] = score.mean()
        dic_community[f_n]['roc_std'] = score.std()
        print("ROC AUC: %.4f (%.4f)" % (score.mean(), score.std()))

    pickle.dump(dic_community, open("svm_results/r_" + react + ".pkl", "wb"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--react_name",
                        help='Reactome pathway on which this code will be executed. For exact names, check meta_data/ReactomeData/ folder')
    args = parser.parse_args()
    print("Going with", args.react_name)

    # Creating the whole table
    all_pandas = []
    for f_name in TISSUES:
        pd_tmp = pd.read_csv("data_filtered/only_geneids_CORRECTED_" + f_name + ".csv", index_col=0)
        pd_tmp.rename(index=lambda x: x + "_" + f_name, inplace=True)

        all_pandas.append(pd_tmp)

    all_df = pd.concat(all_pandas, sort=False)

    # Manually making the scaling because of NaNs
    all_df = all_df.sub(all_df.min()).div((all_df.max() - all_df.min()))

    all_df.fillna(0, inplace=True)


    # Adding tissue information to the dataframe
    def label_race(row):
        splitted = row.name.split('_')
        term = ""
        for i in range(1, len(splitted)):
            term += "_" + splitted[i]
        return term


    all_df['tissue'] = all_df.apply(lambda row: label_race(row), axis=1)

    calculate_react_svm(args.react_name)
