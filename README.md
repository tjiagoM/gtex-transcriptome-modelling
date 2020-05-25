[![DOI](https://zenodo.org/badge/265633374.svg)](https://zenodo.org/badge/latestdoi/265633374)


# Multilayer modelling and analysis of the human transcriptome

*Tiago Azevedo, Giovanna Maria Dimitri, Pietro Lio, Eric R. Gamazon*

This repository contains all the code necessary to run and further extend the experiments presented in the following bioRxiv preprint: [https://doi.org/10.1101/2020.05.21.109082](https://doi.org/10.1101/2020.05.21.109082)

## Abstract

In the present work, we performed a comprehensive intra-tissue and inter-tissue network analysis of the human transcriptome.
We generated an atlas of communities in co-expression networks in each of 49 tissues and evaluated their tissue specificity.
UMAP embeddings of gene expression from the identified communities recovered biologically meaningful tissue clusters, based on tissue organ membership or known shared function.
We developed an approach to quantify the conservation of global structure and estimate the sampling distribution of the distance between tissue clusters via bootstrapped manifolds.
We found not only preserved local structure among clearly related tissues (e.g., the 13 brain regions) but also a strong correlation between the clustering of these related tissues relative to the remaining ones.
Interestingly, brain tissues showed significantly higher variability in community size than non-brain (<img src="https://render.githubusercontent.com/render/math?math=p = 1.55 \times 10^{-4}">).
We identified communities that capture some of our current knowledge about biological processes, but most are likely to encode novel and previously inaccessible functional information.
For example, we found a 17-member community present across all of the brain regions, which shows significant enrichment for the nonsense-mediated decay pathway (adjusted <img src="https://render.githubusercontent.com/render/math?math=p = 1.01 \times 10^{-37}">).
We also constructed multiplex architectures to gain insights into tissue-to-tissue mechanisms for regulation of communities in the transcriptome, including communities that are likely to play a functional role throughout the central nervous system (CNS) and communities that may participate in the interaction between the CNS and the enteric nervous system.
Notably, new gene expression data can be embedded into our models to accelerate discoveries in high-dimensional molecular datasets.
Our study provides a rich resource of co-expression networks, communities, multiplex architectures, and enriched pathways in a broad collection of tissues, to catalyse research into inter-tissue regulatory mechanisms and enable insights into their downstream phenotypic consequences.


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
