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

```r
install.packages("devtools")
library(devtools)
install_github("rezakj/scSeqR")
```

- Download and unzip a publicly available sample [PBMC](https://en.wikipedia.org/wiki/Peripheral_blood_mononuclear_cell) scRNA-Seq data.

```r
setwd("/your/download/directory")
sample.file.url = "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
download.file(url = sample.file.url, 
     destfile = "pbmc3k_filtered_gene_bc_matrices.tar.gz", 
     method = "auto")     
untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")    
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
my.obj <- norm.data(my.obj, 
     "ranked.glsf",
     top.rank = 500) # best for scRNA-Seq

# more examples
#my.obj <- norm.data(my.obj, "global.glsf") # best for bulk RNA-Seq 
#my.obj <- norm.data(my.obj, "rpm", rpm.factor = 1000) # best for bulk RNA-Seq
#my.obj <- norm.data(my.obj, "spike.in", spike.in.factors = NULL)
#my.obj <- norm.data(my.obj, "no.norm") # if the data is already normalized
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
     clust.dim = 2, 
     clust.type = "tsne")

# more examples
#my.obj <- cluster.data(my.obj, clust.method = "base.mean.rank", top.rank = 500, clust.dim = 3, clust.type = "tsne")
#my.obj <- cluster.data(my.obj, clust.method = "base.mean.rank", top.rank = 500, clust.dim = 2, clust.type = "pca")
#my.obj <- cluster.data(my.obj, clust.method = "base.mean.rank", top.rank = 500, clust.dim = 3, clust.type = "pca")
#my.obj <- cluster.data(my.obj, clust.method = "base.mean.rank", top.rank = 500, clust.type = "distance") # nor recomanded for scRNA-Seq
```        
- plot data before cluster assignment.

```r
cluster.plot(my.obj,
     cell.size = 1.75, 
     plot.type = "tsne", 
     clust.assigned = FALSE, 
     clust.dim = 2, 
     cell.color = "blue")

# more examples 
#cluster.plot(my.obj, plot.type = "tsne", clust.assigned = FALSE, clust.dim = 3)
#cluster.plot(my.obj,cell.size = 1.75, plot.type = "pca", clust.assigned = FALSE, clust.dim = 2, cell.color = "blue")
#cluster.plot(my.obj, plot.type = "pca", clust.assigned = FALSE, clust.dim = 3)

# for interactive plots
htmlwidgets::saveWidget(ggplotly(cluster.plot(my.obj, 
     plot.type = "tsne", 
     clust.assigned = FALSE, 
     clust.dim = 3)), "tSNE_plot3d.html")
```

- Find optimal number of clusters          

```r
opt.clust.num(my.obj, max.clust = 10, gap.stat.nboot = 50, clust.type = "tsne", clust.dim = 2, opt.method = "silhouette")
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/optim_clust_num1.png" width="800" />
</p>
              
- Assign clusters

```r
my.obj <- assign.clust(my.obj, 
     clust.num = 7,
     clust.type = "tsne",
     clust.dim = 2)

# more examples
#my.obj <- assign.clust(my.obj, clust.num = 7,clust.type = "tsne",clust.dim = 3)
#my.obj <- assign.clust(my.obj, clust.num = 2,clust.type = "pca",clust.dim = 2)
#my.obj <- assign.clust(my.obj, clust.num = 2,clust.type = "pca",clust.dim = 3)
```

- Plot clusters

```r
cluster.plot(my.obj,
     cell.size = 1.75, 
     plot.type = "tsne",
     clust.assigned = TRUE, 
     clust.dim = 2)

# more examples
#cluster.plot(my.obj,cell.size = 1.75, plot.type = "tsne",clust.assigned = TRUE, clust.dim = 3)
#cluster.plot(my.obj,cell.size = 1.75, plot.type = "pca",clust.assigned = TRUE, clust.dim = 2)
#cluster.plot(my.obj,cell.size = 1.75, plot.type = "pca",clust.assigned = TRUE, clust.dim = 3)

# for interactive plots 
htmlwidgets::saveWidget(ggplotly(cluster.plot(my.obj,
     cell.size = 1.75, 
     plot.type = "tsne",
     clust.assigned = TRUE, 
     clust.dim = 3)), "tSNE_plot3d_clustered.html")
     
     # or
     
htmlwidgets::saveWidget(ggplotly(cluster.plot(my.obj,
     cell.size = 1.75, 
     plot.type = "tsne",
     clust.assigned = TRUE, 
     clust.dim = 2)), "tSNE_plot2d_clustered.html")     
```

To see the above made interactive plots click on these links: [2Dplot](https://genome.med.nyu.edu/results/external/scSeqR/doc/tSNE_plot2d_clustered.html) and [3Dplot](https://genome.med.nyu.edu/results/external/scSeqR/doc/tSNE_plot3d_clustered.html)
        
<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/tSNE_plot.png" width="800" height="700" />
</p>
        

- Save your object

```r
save(my.obj, file = "my.obj.Robj")
```        
     

- Plot genes

```r
gene.plot(my.obj, gene = "MS4A1", box.to.test = 2, box.pval = "sig.signs")
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/master/doc/plot_MS4A1_Pval_signs.png" width="700" height="700" />
  
</p>









