import argparse
import pickle

import numpy as np
import pandas as pd
from sklearn import svm
from sklearn.model_selection import cross_validate

from definitions import TISSUES


# python -u 05_01_svms_communities.py --tissue_num NUM | tee outputs/output_05_01_NUM.txt


def calculate_svm(tissue_name):
    corr_mat = pd.read_pickle("data/corr_" + tissue_name + ".pkl")
    communities = pickle.load(open("results/louvain_modules_" + tissue_name + ".pkl", "rb"))

    community_id = 0

    # For each community found in each tissue
    for community in np.unique(communities[0]):
        common = [x for x in all_df.columns.values if x in np.array(corr_mat.columns)[communities[0] == community]]

        # Communitites of size 3 or less are considered small
        if len(common) <= 3:
            continue

        filtered_df = all_df.loc[:, common]

        filtered_w_tissue = filtered_df.join(all_df['tissue'])

        dic_community = {}
        dic_community['genes'] = common
        community_id += 1

        X = filtered_w_tissue.loc[:, common].values

        for f_name in sorted(TISSUES):
            print("-- Predicting " + f_name + " with " + tissue_name + " community " + str(community_id))
            y = filtered_w_tissue.loc[:, 'tissue'].copy().values
            for j, elem in enumerate(y):
                if elem == "_" + f_name:
                    y[j] = 1
                else:
                    y[j] = 0

            dic_community[f_name] = {}

            scoring = ['accuracy', 'f1', 'roc_auc']

            clf = svm.SVC(kernel='linear', class_weight="balanced")
            scores = cross_validate(clf, X, list(y), cv=3, scoring=scoring, n_jobs=-1)

            score = scores['test_accuracy']
            dic_community[f_name]['acc'] = score.mean()
            dic_community[f_name]['acc_std'] = score.std()
            print("Accuracy: %.4f (%.4f)" % (score.mean(), score.std()))

            score = scores['test_f1']
            dic_community[f_name]['f1'] = score.mean()
            dic_community[f_name]['f1_std'] = score.std()
            print("F1 score: %.4f (%.4f)" % (score.mean(), score.std()))

            score = scores['test_roc_auc']
            dic_community[f_name]['roc'] = score.mean()
            dic_community[f_name]['roc_std'] = score.std()
            print("ROC AUC: %.4f (%.4f)" % (score.mean(), score.std()))

        pickle.dump(dic_community, open("svm_results/" + tissue_name + '_' + str(community_id) + ".pkl", "wb"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--tissue_num", type=int)
    args = parser.parse_args()
    print("Going with", TISSUES[args.tissue_num])

    # Creating the whole table
    all_pandas = []
    for f_name in sorted(TISSUES):
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

    calculate_svm(TISSUES[args.tissue_num])
