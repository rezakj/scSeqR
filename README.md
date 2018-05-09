# scSeqR

### Single Cell Sequencing R package (scSeqR)

scSeqR is an R package that can analyze single sequencing data types (i.e scRNA-Seq) and large matrix files. The program inputs single cell data in 10X format or standard matrix and data frames and helps you to perform QC, filtering, visualization, normalization, clustering, differential expression analysis and find positive and negative markers for each cluster.

### How to install scSeqR

- First clone the package.

        git clone https://github.com/rezakj/scSeqR.git

- Install the dependencies for scSeqR.

        install.packages(c("ggplot2",
        "Matrix",
        "Rtsne",
        "gmp", 
        "factoextra", 
        "gridExtra"0))
        
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


        plot(my.obj@stats$UMIs,
        my.obj@stats$mito.percent, 
        main = "UMIs/Mito", 
        ylab = "Percent Mito", 
        xlab = "UMIs")
        
        plot(my.obj@stats$UMIs,
        my.obj@stats$nGenes, 
        main = "UMIs/genes", 
        ylab = "genes", 
        xlab = "UMIs")

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/UMIs_Mito.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/UMIs_genes.png" width="400"/>      
</p>


- Filter cells

        my.obj <- filter.cells(my.obj,
                min.mito = 0, 
                max.mito = 0.05, 
                min.genes = 200, 
                max.genes = 2500, 
                min.umis = 0, 
                max.umis = Inf)
                
        dim(my.obj@main.data)


- Normalize data 

        my.obj <- norm(my.obj, "ranked.glsf", top.rank = 500)

- Scale data 

        my.obj <- scale.data(my.obj)

- Cluster data 

        my.obj <- cluster.data(my.obj, 
                clust.method = "base.mean.rank", 
                top.rank = 500, 
                clust.dim = 2)
        
- Find optimal number of clusters          

        pot.clust.num(my.obj)
        
<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/optim_clust_num1.png" width="800"/>
</p>
              
- Assign clusters and plot them

        my.obj <- assign.clust(my.obj, clust.num = 7)
        tsne.plot(my.obj)

        
<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/tSNE_plot.png" width="800"/>
</p>
        

- Make intractive plot in html format

        library(plotly)
        htmlwidgets::saveWidget(ggplotly(tsne.plot(my.obj)), "tSNE_plot.html")


- Save your object
        save(my.obj, file = "my.obj.Robj")
        
     

- Plot genes

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/tSNE_CD14.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/boxplot_CD14.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/tSNE_CD3E.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/boxplot_CD3E.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/tSNE_FCGR3A.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/boxplot_FCGR3A.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/tSNE_GNLY.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/boxplot_GNLY.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/tSNE_MS4A1.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/boxplot_MS4A1.png" width="400"/>
  
</p>









