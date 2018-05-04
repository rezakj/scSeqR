# scSeqR

### Single Cell RNA-Seq R package (scSeqR)

scSeqR is an R package that can analyze single cell RNA-Seq and large matrix files. The program inputs single cell data in 10X format or standard matrix and data frames and helps you to perform QC, filtering, visualization, normalization, clustering, differential expression analysis and find positive and negative markers for each cluster.

### How to install scSeqR

- First clone the package.

        git clone https://github.com/rezakj/scSeqR.git

- Install the dependencies for scSeqR.

        install.packages(c("ggplot2","Matrix","Rtsne"))

- Then install the package in R.

        install.packages('directory/to/scSeqR', repos = NULL, type="source")

- Then load the package in R.

        library("scSeqR")


### How to use scSeqR

To run a test sample follow these steps:

- Download this data from here (using terminal).

        wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

- Unzip the data.

        tar xvzf pbmc3k_filtered_gene_bc_matrices.tar.gz

- Go to the R environment load the scSeqR package and the PBMC data.

        library("scSeqR")
        my.data <- load10x("filtered_gene_bc_matrices/hg19/",gene.name = "geneSymbol")

To see the help page for each function use question mark as: 

        ?load10x

There is also a sample data that comes with the package. To see the head and the structure of the sample data issue this command:

       head(sample.10x.data)[1:5]

You can load the sample data as:

        my.data <- sample.10x.data
        
        
- Make an object of class scSeqR

        my.obj <- make.obj(my.data)

- Perform some QC 

        dim(my.obj@raw.data)
        my.obj <- UMIs.genes.mito(my.obj)
        summary(my.obj@stats)
        stats.plot(my.obj)


<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/stats.png" width="800"/>
</p>
