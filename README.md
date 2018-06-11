# scSeqR

### official release in July (under construction)

### Single Cell Sequencing R package (scSeqR)

scSeqR is an R package that can analyze single cell sequencing data types (i.e [scRNA-seq](https://en.wikipedia.org/wiki/Single_cell_sequencing#Single-cell_RNA_sequencing_(scRNA-seq))) and large numeric [matrix](https://en.wikipedia.org/wiki/Matrix_(mathematics)) 
files (i.e. count tables with many samples from [TCGA](https://cancergenome.nih.gov/)). The program inputs single cell data in [10X format](https://www.10xgenomics.com/), large numeric **matrix files** and **data frames** and helps you to perform QC, filtering, visualization, normalization, clustering, differential expression analysis and find positive and negative markers for each cluster. scSeqR, allows you to choose from **multiple normalization** methods and **spike-in normalization** depending on your data type. Alternatively, you can also use 
**already normalized** data.

***
## How to install scSeqR

- Install the dependencies for scSeqR in R.

```r
install.packages(c("ggplot2",
     "Matrix",
     "plotly",
     "Rtsne",
     "gmp", 
     "factoextra", 
     "gridExtra"))
 ```
        
- Then install the package in R.

```r
library(devtools)
install_github("rezakj/scSeqR")
```

## Download a sample data

- Download and unzip a publicly available sample [PBMC](https://en.wikipedia.org/wiki/Peripheral_blood_mononuclear_cell) scRNA-Seq data.

```r
# set your working directory 
setwd("/your/download/directory")

# save the URL as an object
sample.file.url = "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"

# download the file
download.file(url = sample.file.url, 
     destfile = "pbmc3k_filtered_gene_bc_matrices.tar.gz", 
     method = "auto")  

# unzip the file. 
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

- Aggregating data
     
Conditions in scSeqR, are set in the header of the data and are separated by an underscore (_).
Let's say you want to merge multiple datasets and run scSeqR in aggregated mode. To do this let's divide your sample into 3,  assuming that there are 3 samples and aggregate them into one matrix. 

```r
dim(my.data)
# [1] 32738  2700

# divide your sample into three samples for this example 
  sample1 <- my.data[1:900]
  sample2 <- my.data[901:1800]
  sample3 <- my.data[1801:2700]
  
# merge all of your samples to make a single aggregated file.    
my.data <- data.aggregation(samples = c("sample1","sample2","sample3"), 
	condition.names = c("WT","KO","Ctrl"))
```

- Check the head of your file.

```r
# here is how the head of the first 2 cells in the aggregated file looks like.	
head(my.data)[1:2]
#         WT_AAACATACAACCAC-1 WT_AAACATTGAGCTAC-1
#A1BG                       0                   0
#A1BG.AS1                   0                   0
#A1CF                       0                   0
#A2M                        0                   0
#A2M.AS1                    0                   0

# as you see the header has the conditions now
```


- Make an object of class scSeqR.

```r
my.obj <- make.obj(my.data)
```

- Perform some QC 

```r
my.obj <- UMIs.genes.mito(my.obj)
dim(my.obj@raw.data)
# [1] 32738  2700
``` 

- Plot QC

By defult all the plotting functions would creat interactive html files unless interactive = FALSE.

```r
# plot UMIs, genes and percent mito all at once and in one plot. 
# you can make them individually as well, see the arguments ?stats.plot.
stats.plot(my.obj,
	plot.type = "box.gene.umi.mito",
	out.name = "UMI-plot",
	interactive = FALSE,
	cell.color = "slategray3", 
	cell.size = 1, 
	cell.transparency = 0.5,
	box.color = "red",
	box.line.col = "green")
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/stats.png" width="800" height="700" />
</p>

```r  
# Scatter plots
stats.plot(my.obj, plot.type = "point.mito.umi", out.name = "mito-umi-plot")
stats.plot(my.obj, plot.type = "point.gene.umi", out.name = "gene-umi-plot")
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/mito.umi.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/gene.umi.png" width="400"/>      
</p>

To see an example interactive plot click on this links: [mito-UMIs plot](lllll)

- Filter cells. 

scSeqR allows you to filter based on library sizes (UMIs), number of genes per cell, percent mitochondial content and even based on one or more genes. 

Befor filtering let's do a quick test and see how the data looks like. This might be very helpful to check the gene or genes of your interest and see in howmany cells they expressed. This could also be helful to choose the genes for cell sorting or filteering that are needed in some studies. 

```r
# run gene stats 
my.obj <- gene.stats(my.obj, which.data = "raw.data")

# sort and see the head of the gene stats
head(my.obj@gene.data[order(my.obj@gene.data$numberOfCells, decreasing = T),])
#       genes numberOfCells totalNumberOfCells percentOfCells  meanExp      SDs
#30303 TMSB4X          2700               2700      100.00000 46.00370 30.86793
#3633     B2M          2699               2700       99.96296 44.94926 25.89472
#14403 MALAT1          2699               2700       99.96296 59.88333 32.56044
#27191 RPL13A          2698               2700       99.92593 28.46037 16.70928
#27185  RPL10          2695               2700       99.81481 32.78407 17.90696
#27190  RPL13          2693               2700       99.74074 28.55963 17.00709
```

Let's say that you are only interested in the cells in which RPL13 "OR" RPL10 are expressed, while have a maximum mito percent of 0.05 and a minimum number of 250 genes expressed and a maximum of 2400 genes. 

```r
my.obj <- cell.filter(my.obj,
	filter.by.gene = c("RPL13","RPL10"),
	min.mito = 0,
	max.mito = 0.05,
	min.genes = 200,
	max.genes = 2400,
	min.umis = 0,
	max.umis = Inf)

# chack to see how many cells are left.  
dim(my.obj@main.data)
# [1] 32738  2634
```

- Normalize data

You have a few options to normalize your data based on your study. You can also normalize your data using any tool other than scSeqR. We recomend "ranked.glsf" normalization for most singel cell studies, this is Geometric Library Size Factor (GLSF) normalization that is using top expressed genes ranked by base mean. This normalization is great for fixing for matrices with a lot of zeros and because it's geometric it is great for fixing for batch effects as long as all the data is aggregated in one file (to aggregate your data see "aggregating data" section above). 

```r
my.obj <- norm.data(my.obj, 
     "ranked.glsf",
     top.rank = 500) # best for scRNA-Seq

# more examples
#my.obj <- norm.data(my.obj, "global.glsf") # best for bulk RNA-Seq 
#my.obj <- norm.data(my.obj, "rpm", rpm.factor = 100000) # best for bulk RNA-Seq
#my.obj <- norm.data(my.obj, "spike.in", spike.in.factors = NULL)
#my.obj <- norm.data(my.obj, "no.norm") # if the data is already normalized
```

- Scale data 

This would greatly help in clustring and plotting. 

```r
my.obj <- scale.data(my.obj)
```
- Gene stats

It's better to run gene.stats on your main data and update the gene.data of your object. 

```r
my.obj <- gene.stats(my.obj, which.data = "main.data")

head(my.obj@gene.data[order(my.obj@gene.data$numberOfCells, decreasing = T),])
#       genes numberOfCells totalNumberOfCells percentOfCells  meanExp      SDs
#30303 TMSB4X          2634               2634      100.00000 46.74179 24.85859
#3633     B2M          2633               2634       99.96203 46.98768 25.58558
#14403 MALAT1          2633               2634       99.96203 65.53506 39.68695
#27191 RPL13A          2633               2634       99.96203 28.96777 12.96873
#27185  RPL10          2632               2634       99.92407 32.74179 11.13561
#27190  RPL13          2630               2634       99.84814 29.12121 13.73905
```

- Make a gene model for clustering

It's best to always to avoid global clustering and use a set of model genes. In bulk RNA-seq data it is very common to cluster the samples based on top 500 genes ranked by base mean, this is to reduce the noise. In scRNA-seq data, it's great to do so as well. This coupled with our ranked.glsf normalization is great for matrices with a lot of zeros. You can also use your set of genes as a model rather than making one. 

```r
make.gene.model(my.obj, 
	dispersion.limit = 1.5, 
	base.mean.rank = 500, 
	no.mito.model = T, 
	mark.mito = T, 
	interactive = T,
	out.name = "gene.model")
```
To view an the html intractive plot click on this links: [Dispersion plot](https://rawgit.com/rezakj/scSeqR/dev/doc/gene.model.html)


<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/gene.model.png" width="800" height="800" />
</p>


- Cluster data

scSeqR let's you perform 2 dimentional and 3 dementional PCA and tSNE clustering and you also have the option of distance based clustering using "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski" methods. 

It's best to use tSNE for singel cell data and large matrices.  

```r
# for 2D tSNE
my.obj <- cluster.data(my.obj, 
	clust.method = "gene.model", 
	gene.list = "my_model_genes.txt", 
	clust.dim = 2, 
	clust.type = "tsne")
# for 2D PCA	
my.obj <- cluster.data(my.obj, 
	clust.method = "gene.model", 
	gene.list = "my_model_genes.txt", 
	clust.dim = 2, 
	clust.type = "pca")	

# more examples
#my.obj <- cluster.data(my.obj, clust.method = "base.mean.rank", top.rank = 500, clust.dim = 3, clust.type = "tsne")
#my.obj <- cluster.data(my.obj, clust.method = "base.mean.rank", top.rank = 500, clust.dim = 2, clust.type = "pca")
#my.obj <- cluster.data(my.obj, clust.method = "base.mean.rank", top.rank = 500, clust.dim = 3, clust.type = "pca")
#my.obj <- cluster.data(my.obj, clust.method = "base.mean.rank", top.rank = 500, clust.type = "distance") # nor recomanded for scRNA-Seq
```        

- Visualize conditions

```r
# tSNE
cluster.plot(my.obj,
	plot.type = "tsne",
	col.by = "conditions",
	clust.dim = 2,
	interactive = F)
# pca 
cluster.plot(my.obj,
	plot.type = "tsne",
	col.by = "conditions",
	clust.dim = 2,
	interactive = F)
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/tSNE_conditions.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/PCA_conditions.png" width="400"/>      
</p>

- Find optimal number of clusters

scSeqR allows to choose from 3 different algorisms ("elbow.wss", "silhouette", "gap.stat") to find the optimal number of clusters. 

```r
opt.clust.num(my.obj, max.clust = 20, 
	clust.type = "tsne", 
	clust.dim = 2, 
	opt.method = "elbow.wss")

# more examples 
#opt.clust.num(my.obj, max.clust = 20, gap.stat.nboot = 200, clust.type = "tsne", clust.dim = 2, opt.method = "gap.stat")
#opt.clust.num(my.obj, max.clust = 10, clust.type = "tsne", clust.dim = 2, opt.method = "silhouette")
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/optim_clust_num1.png" width="800" />
</p>
   
- Assign clusters

```r
my.obj <- assign.clust(my.obj, 
     clust.num = 4,
     clust.type = "tsne",
     clust.dim = 2)
```

- Plot clustered data 

```r
cluster.plot(my.obj,
	cell.size = 1,
	plot.type = "tsne",
	cell.color = "black",
	back.col = "white",
	col.by = "clusters",
	cell.transparency = 0.5,
	clust.dim = 2,
	interactive = T,
	out.name = "tSNE")

# more examples
#cluster.plot(my.obj,cell.size = 1.75, plot.type = "tsne",clust.assigned = TRUE, clust.dim = 3)
#cluster.plot(my.obj,cell.size = 1.75, plot.type = "pca",clust.assigned = TRUE, clust.dim = 2)
#cluster.plot(my.obj,cell.size = 1.75, plot.type = "pca",clust.assigned = TRUE, clust.dim = 3)   
```

To see the above made interactive plots click on these links: [2Dplot](https://rawgit.com/rezakj/scSeqR/master/doc/tSNE_plot2d_clustered.html) and [3Dplot](https://rawgit.com/rezakj/scSeqR/master/doc/tSNE_plot3d_clustered.html)
        
<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/tSNE.png" width="700" height="700" />
</p>
        
- Avrage expression per cluster

```r
my.obj <- clust.avg.exp(my.obj, clust.type = "tsne", clust.dim = 2)
head(my.obj@clust.avg)
```

|gene	|cluster_1	|cluster_2	|cluster_3	|cluster_4	|cluster_5	|cluster_6	|cluster_7|
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
|A1BG	|0.0916392944	|0.0436354813	|0.0343642627	|0.0667622253	|0.0613617813	|0.0746466649	|0.0840048157|
|A1BG.AS1	|0.0110216269	|0.0149905995	|0.0121888084	|0.0023281341	|0.0076570435	|0.0151865099	|0.0099459357|
|A1CF	|0	|0	|0	|0	|0	|0	|0|
|A2M	|0	|0	|0.0027540305	|0.0023446032	|0.0068563104	|0.0013456149	|0|
|A2M.AS1	|0.0204136324	|0	|0.0022296427	|0.04610626	|0.0128255855	|0.0064323158	|0.0054419049|


- Save your object

```r
save(my.obj, file = "my.obj.Robj")
```        
     

- Plot genes

```r
gene.plot(my.obj, gene = "MS4A1", box.to.test = 2, box.pval = "sig.signs")
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/plot_MS4A1_Pval_signs.png" width="700" height="700" />
  
</p>

