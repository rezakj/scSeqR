# scSeqR

### Single Cell Sequencing R package (scSeqR)

scSeqR is an R package that can analyze single sequencing data types (i.e [scRNA-seq](https://en.wikipedia.org/wiki/Single_cell_sequencing#Single-cell_RNA_sequencing_(scRNA-seq))) and large numeric [matrix](https://en.wikipedia.org/wiki/Matrix_(mathematics)) 
files (i.e. count tables with many samples from [TCGA](https://cancergenome.nih.gov/)). The program inputs single cell data in [10X format](https://www.10xgenomics.com/), large numeric **matrix files** and **data frames** and helps you to perform QC, filtering, visualization, normalization, clustering, differential expression analysis and find positive and negative markers for each cluster. scSeqR, allows you to choose from **multiple normalization** methods and **spike-in normalization** depending on your data type. Alternatively, you can also use 
**already normalized** data.

***
## How to install scSeqR

- Install the dependencies for scSeqR in R.

```r
install.packages(c("ggplot2",
     "Matrix",
     "Rtsne",
     "gmp", 
     "factoextra", 
     "gridExtra"))
 ```
        
- Then install the package in R.

```console
install.packages("devtools")
library(devtools)
install_github("rezakj/scSeqR")
```

- Download and unzip a publicly available sample [PBMC](https://en.wikipedia.org/wiki/Peripheral_blood_mononuclear_cell) scRNA-Seq data.

```console
# shell (bash) 
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

tar xvzf pbmc3k_filtered_gene_bc_matrices.tar.gz
```

***
# How to use scSeqR for analyzing scRNA-seq data

To run a test sample follow these steps:

- Go to the R environment load the scSeqR package and the PBMC sample data that you downloaded.

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
head(sample.10x.data)[1:4]
```

|  | AAACATACAACCAC | AAACATTGAGCTAC | AAACATTGATCAGC | AAACCGTGCTTCCG |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| MALAT1        |     49        |    142        |    171        |     11        |
| TMSB4X        |     47        |     62        |    117        |    114        |
| B2M           |     76        |     75        |     69        |     41        |
| RPL10         |     34        |     92        |     49        |     22        |
| RPL13         |     29        |     45        |     16        |     15        |
| RPL13A        |     37        |     81        |     40        |     16        |

        
Conditions in scSeqR are set in the header of the data and are separated by an underscore (_) as below:

|  | condition1_AAACATACAACCAC | condition1_AAACATTGAGCTAC | ... | condition2_AAACATTGATCTGC | condition2_AAACCGTGCTTGCG |
| ------------- | ------------- | ------------- | ------------- | ------------- |------------- |
| MALAT1        |     49        |    142        |    ...        |     112       |     100      |
| TMSB4X        |     47        |     62        |    ...        |    11         |     81       |
| B2M           |     76        |     75        |    ...        |     45        |     51       |
| RPL10         |     34        |     92        |    ...        |     26        |     18       |
| RPL13         |     29        |     45        |    ...        |     75        |     110      |
| RPL13A        |     37        |     81        |    ...        |     66        |     12       |

Follow [this] tutorial to learn how to merge multiple datasets and run scSeqR in aggregated mode.

- Make an object of class scSeqR

```r
my.obj <- make.obj(my.data)
```

- Perform some QC 

```r
dim(my.obj@raw.data)
```

| # of genes/rows| # of cells/columns |
| ------------- | ------------- |
| 32738        |     2700        |

```r
my.obj <- UMIs.genes.mito(my.obj)
stats.plot(my.obj)
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/stats.png" width="800" height="700" />
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

| # of genes/rows| # of cells/columns |
| ------------- | ------------- |
| 32738        |     2638        |

- Normalize data 

```r
my.obj <- norm(my.obj, norm.method = "ranked.glsf", top.rank = 500)
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
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/optim_clust_num1.png" width="800" />
</p>
              
- Assign clusters and plot them

```r
my.obj <- assign.clust(my.obj, clust.num = 7)
tsne.plot(my.obj)
```
        
<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/tSNE_plot.png" width="800" height="700" />
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
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/plot_CD14.png" width="800" height="700" />
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/plot_MS4A1.png" width="800" height="700" />
  
</p>









