[![DOI](https://zenodo.org/badge/265633374.svg)](https://zenodo.org/badge/latestdoi/265633374)


# Multilayer modelling and analysis of the human transcriptome

*Tiago Azevedo, Giovanna Maria Dimitri, Pietro Lio, Eric R. Gamazon*

This repository contains all the code necessary to run and further extend the experiments presented in the following bioRxiv preprint: [https://doi.org/10.1101/2020.05.21.109082](https://doi.org/10.1101/2020.05.21.109082)

## Abstract

Here, we performed a comprehensive intra-tissue and inter-tissue network analysis of the human transcriptome.
We generated an atlas of communities in co-expression networks in 49 tissues (GTEx v8), evaluated their tissue specificity, and investigated their methodological implications.
UMAP embeddings of gene expression from the communities (representing nearly 18% of all genes) robustly identified biologically-meaningful clusters.
Methodologically, integration of the communities into a transcriptome-wide association study of C-reactive protein (CRP) in 361,194 individuals in the UK Biobank identified genetically-determined expression changes associated with CRP and led to considerably improved performance.
Furthermore, a deep learning framework applied to the communities in nearly 11,000 tumours profiled by The Cancer Genome Atlas across 33 different cancer types learned biologically-meaningful latent spaces, representing metastasis (<img src="https://render.githubusercontent.com/render/math?math=p < 2.2 \times 10^{-16}">) and stemness (<img src="https://render.githubusercontent.com/render/math?math=p < 2.2 \times 10^{-16}">).
Our study provides a rich genomic resource to catalyse research into inter-tissue regulatory mechanisms and their downstream phenotypic consequences.


## Repository Structure

The number of each script in the root of this repository corresponds to the order in which they are run in the paper.


Each folder used in this repository is explained as follows:

* `meta_data`: This folder includes some completementary files to GTEx necessary to run some experiments. Examples of such files include phenotype information, as well as information about conversion of gene names and identification of reactome pathways.

* `outputs`: This folder contains the outputs of some of the numbered scripts. Some contain important information used in the paper, as it is, for example, the files `output_02_01.txt`, `output_04_02.txt`, and `output_06_02.txt`.

* `results`: Some result files, like the communities identified from the Louvain algorithm.

* `svm_results`: Files with the metrics resulting from the SVM predictions.


In the repository there are also some jupyter notebooks which we hope can help researchers in using our results in their own experiments, as well as improve the reproducibility of this paper:

* `09_community_info.ipyb`: Instructions on how to check information regarding each community characterisation, including generation of LaTeX code.

* `10_reactomes_per_tissue.ipynb`: Instructions on how to check information regarding which reactomes were able to predict each tissue.

* `11_multiplex_enrichment.ipynb`: Instructions on how to check the group of genes identified in each multiplex network.

* `12_tcga.ipynb`: the code used in the paper to analyse the TCGA dataset within the GTEx pipeline of the paper, as well as a targeted R code (`12_01_correct_confounds_tcga.R`) used to correct the data.

* `13_plots_for_paper.ipynb`: the code used to generate the plots from the paper.
