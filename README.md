# scSeqR

### Single Cell Sequencing R package (scSeqR)

scSeqR is an R package that can analyze single sequencing data types (i.e [scRNA-seq](https://en.wikipedia.org/wiki/Single_cell_sequencing#Single-cell_RNA_sequencing_(scRNA-seq)) :+1:) and large [matrix](https://en.wikipedia.org/wiki/Matrix_(mathematics)) :+1:
files (i.e. count tables with many samples from [TCGA](https://cancergenome.nih.gov/)). The program inputs single cell data in [10X format](https://www.10xgenomics.com/) :+1: or standard matrix and data frames and helps you to perform QC, filtering, visualization, normalization, clustering, differential expression analysis and find positive and negative markers for each cluster. scSeqR, allows you to choose from multiple normalization methods :+1: and spike-in normalization :+1: depending on your data type. Alternatively, you can also use already normalized :+1: data.

***
## How to install scSeqR

- First clone the package.

```console
# shell (bash)
git clone https://github.com/rezakj/scSeqR.git
```
- Install the dependencies for scSeqR.

```r
install.packages(c("ggplot2",
     "Matrix",
     "Rtsne",
     "gmp", 
     "factoextra", 
     "gridExtra"))
 ```
        
- Then install the package in R.

```r
install.packages('directory/to/scSeqR', repos = NULL, type="source")
```
- Then load the package in R.

```r
library("scSeqR")
```

***
# How to use scSeqR

To run a test sample follow these steps:

- Download sample data and unzip.

```console
# shell (bash) 
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

tar xvzf pbmc3k_filtered_gene_bc_matrices.tar.gz
```

- Go to the R environment load the scSeqR package and the PBMC data.

```r
library("scSeqR")
my.data <- load10x("filtered_gene_bc_matrices/hg19/",gene.name = "geneSymbol")
```

To see the help page for each function use question mark as: 

```r
?load10x
```

There is also a sample data that comes with the package. To see the head and the structure of the sample data issue this command:

```r
head(sample.10x.data)[1:5]
```

| First Header  | Second Header |
| ------------- | ------------- |
| Content Cell  | Content Cell  |
| Content Cell  | Content Cell  |
  
        
- Make an object of class scSeqR

```r
my.obj <- make.obj(my.data)
```

- Perform some QC 

```r
dim(my.obj@raw.data)
my.obj <- UMIs.genes.mito(my.obj)
summary(my.obj@stats)
stats.plot(my.obj)
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/stats.png" width="800"/>
</p>

```r
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
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/UMIs_Mito.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/UMIs_genes.png" width="400"/>      
</p>


- Filter cells

```r
my.obj <- filter.cells(my.obj,
     min.mito = 0, 
     max.mito = 0.05, 
     min.genes = 200, 
     max.genes = 2500, 
     min.umis = 0, 
     max.umis = Inf)
                
dim(my.obj@main.data)
```

- Normalize data 

```r
my.obj <- norm(my.obj, "ranked.glsf", top.rank = 500)
```

- Scale data 

```r
my.obj <- scale.data(my.obj)
```

- Cluster data

```r
my.obj <- cluster.data(my.obj, 
       clust.method = "base.mean.rank", 
       top.rank = 500, 
       clust.dim = 2)
```        
- Find optimal number of clusters          

```r
pot.clust.num(my.obj)
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/optim_clust_num1.png" width="800"/>
</p>
              
- Assign clusters and plot them

```r
my.obj <- assign.clust(my.obj, clust.num = 7)
tsne.plot(my.obj)
```
        
<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/tSNE_plot.png" width="800"/>
</p>
        

- Make intractive plot in html format

```r
library(plotly)
htmlwidgets::saveWidget(ggplotly(tsne.plot(my.obj)), "tSNE_plot.html")
```

- Save your object

```r
save(my.obj, file = "my.obj.Robj")
```        
     

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









