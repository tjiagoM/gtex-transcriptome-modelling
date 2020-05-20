import pickle
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable

from definitions import TISSUES

indexes = []
all_data_f1 = []

reactomes = pd.read_table("meta_data/reactomes_names.txt", header=None).iloc[:, 0].values

for react in reactomes:
    dic_community = pickle.load(open("svm_results/r_" + react + ".pkl", "rb"))

    # Going over each tissue
    arr_f1 = []

    for f_name in TISSUES:
        arr_f1.append(dic_community[f_name]['f1'])

    all_data_f1.append(arr_f1)
    indexes.append(react)

all_f1 = pd.DataFrame(all_data_f1, index=indexes, columns=TISSUES)

##########
# Plotting it all
fig, ax = plt.subplots(figsize=(60, 120))

im = ax.imshow(all_f1.values, cmap='bwr', interpolation='none', vmin=0, vmax=1)

ax.set_xticks(np.arange(len(TISSUES)))
ax.set_xticklabels(TISSUES)

ax.set_yticks(np.arange(len(indexes)))
ax.set_yticklabels(indexes)

plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

ax.set_title("F1 scores of each REACTOME predicting tissue", fontweight='bold')
ax.set_ylabel("Reactomes", fontweight='bold')
ax.set_xlabel("Tissues being predicted", fontweight='bold')

# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.5 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.5)

fig.colorbar(im, cax=cax)
fig.tight_layout()
plt.savefig("figures/f1_reactomes_tissues.png")
plt.close(fig)

##############
# Printing/plotting some other details

# enumerate here just because of order in the folder
for ind, tissue in enumerate(TISSUES):
    print("***** " + tissue + " *****")

    tmp = pd.DataFrame(all_f1.loc[all_f1[tissue] > 0.8, tissue].sort_values(ascending=False).round(decimals=4))
    if tmp.shape[0] == 0:
        print("  Not predicted by any reactome!")
        print()
        continue

    print("Predicted by " + str(tmp.shape[0]) + " reactomes")
    print()
    plt.figure()
    tablea = plt.table(cellText=tmp.values, loc='center', rowLabels=tmp.index.values,
                       colLabels=['F1 Score for ' + tissue])
    tablea.properties()['celld'][(0, 0)].set_text_props(fontweight='bold')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig("figures/react_" + str(ind) + "_" + tissue + ".png", transparent=True, bbox_inches='tight')
    plt.close()

print()
print()
print()

results_per_reactome = defaultdict(list)

for ind in range(all_f1.shape[0]):
    row = all_f1.iloc[ind, :]
    results_per_reactome[sum(row > 0.8)].append(row.name)

for no in sorted(results_per_reactome.keys(), reverse=True):
    print("****** Reactomes which can predict " + str(no) + " tissues ******")
    print("\n".join(sorted(results_per_reactome[no])))
    print()
    print()
    print()
    print()
