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

This repository contains all the scripts which were used in the paper. The number in each script's name (in the root of this repository) corresponds to the order in which they are run in the paper.

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


## Installing Dependencies

These scripts were tested in a **Linux Ubuntu 16.04** operating system, with environments created using **Anaconda**. 

### Python scripts

We include a working dependency file in `python_environment.yml` describing the exact dependencies used to run the python scripts. In order to install all the dependencies automatically with Anaconda, one can easily just run the following command in the terminal to create an Anaconda environment:

```bash
$ conda env create --force --file python_environment.yml
$ conda activate gtex-env
```

To summarise the `python_environment.yml` file, the main dependencies needed to run these scripts are:

* **gseapy** 0.9.16
* **jupyterlab** 1.1.4
* **matplotlib** 3.1.0
* **networkx** 2.4
* **numpy** 1.17.3
* **pandas** 1.0.1
* **python** 3.7.5
* **scikit-learn** 0.21.3
* **statsmodels** 0.10.2
* **umap-learn** 0.4.2
* **bctpy** 0.5.0

### R scripts

We used R to run the unsupervised correction package *sva*, which is briefly described in the paper. To make things easier, we also include the R dependencies in which these scripts were run, with Anaconda. Similarly to python, one can install them using the following commands:

```bash
$ conda env create --force --file r_environment.yml
$ conda activate r_env
```

After the R environment is created and activated, one should install the *sva* package as descibred in the [original repository](https://bioconductor.org/packages/release/bioc/html/sva.html):

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("sva")
```

We decided to keep python and R scripts in separate environments to avoid dependency issues given they are distinct programming languages.


### Data Requirements

To details on the data used, which cannot be publicly shared in this repository, please see the paper.

## Running the scripts

In order to analyse and see our jupyter notebooks, one just needs to start the jupyter engine in the root of this repository:

```bash
$ jupyter lab --port=8895
```

This command will print a link in the local machine at port 8895, which can be accessed using a browser.

The python scripts in this repository are numbered, suggesting an order by which they should be executed. However, each python script contains at the beginning of the file a small documentation explaining what it does. To run python scripts from `01` to `04_02`, and `06_02`, one just needs to run the following command:

```bash
$ python -u PYTHON_FILE | tee outputs/output_file.txt
```

The previous command will run `PYTHON_FILE` and log the output of the script in *outputs/output_file.txt*. All the other python scripts expect one or two flags to be passed. Information about each flag can be seen in each `parser.add_argument` command in each file, which contains a small documentation of what it means. For example, python script `05_01` expects the flag `--tissue_num`; therefore, that flag needs to be passed when executing the script:

```bash
$ python -u 05_01_svms_communities.py --tissue_num NUM | tee outputs/output_05_01_NUM.txt
```

where `NUM` corresponds to the value for that flag.


The following command is an example of how to run an R script:

```bash
$ Rscript --no-save --no-restore --verbose 02_01_correct_confounds.R > outputs/output_02_01.txt 2>&1
```

The previous command will run the `02_01_correct_confounds.R` script and log the output of the script in *outputs/output_02_01.txt*.


## Other scripts

The scripts for the multilayer modeling approach to TWAS/PrediXcan (CRP in UKB) and Variational Autoencoder model (TCGA) are in this [external repository](https://github.com/gamazonlab/MultilayerModelingTranscriptome).


## Questions?

If you run into any problem or have a question, please just open an issue.
