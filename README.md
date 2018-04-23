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

- Perform some QC 

        my.qc <- UMIs.genes.mito(my.data)
        
        names(my.qc)
        #[1] "mito.percent" "nGenes"       "UMIs"
        
        mito.percent.plot <- ggplot(as.data.frame(my.qc$mito.percent),aes(y=as.data.frame(my.qc$mito.percent),x="mydat", alpha = 0.5)) + 
        theme_bw() + 
        geom_jitter(color = "red") + 
        geom_boxplot() + xlab("my.data") + ylab("percent of mito genes per cell")
        
        nGenes.plot <- ggplot(as.data.frame(my.qc$nGenes),aes(y=as.data.frame(my.qc$nGenes),x="mydat", alpha = 0.5)) + 
        theme_bw() + 
        geom_jitter(color = "red") + 
        geom_boxplot() + xlab("my.data") + ylab("number of genes per cell")
        
        UMIsplot <- ggplot(as.data.frame(my.qc$UMIs),aes(y=as.data.frame(my.qc$UMIs),x="mydat", alpha = 0.5)) + 
        theme_bw() + 
        geom_jitter(color = "red") + 
        geom_boxplot() + xlab("my.data") + ylab("number of UMIs per cell")
        
        pdf("plot_number_of_UMIs.pdf",width = 4,height = 10)
        UMIsplot
        dev.off()
        
        pdf("plot_number_of_genes.pdf",width = 4,height = 10)
        nGenes.plot
        dev.off()
        
        pdf("plot_percent_mito.pdf",width = 4,height = 10)
        mito.percent.plot
        dev.off()

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/plot_number_of_UMIs.png" width="350"/>
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/plot_number_of_genes.pdf" width="350"/>
</p>
