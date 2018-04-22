# scSeqR

### Single Cell RNA-Seq R package (scSeqR)

scSeqR is an R package to analyze Single Cell RNA-Seq and large matrix files. It inputs single cell data in 10X format or standared matrix and data frames and helps you to perform QC, filter, visualize, normalize, cluster and find positive and negetive markers for each cluser using differential expression analysis. 

### How to install scSeqR

- Firts clone the data.

`git clone https://github.com/rezakj/scSeqR.git`

- Then install the package in R.

`install.packages('~/scSeqR', repos = NULL, type="source")`

- Then load the package in R.

`library("scSeqR")`


### How to use scSeqR

To run a test sample follow these steps:

- Download this data from here (using terminal).

`wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz`

- Unzip the data.
