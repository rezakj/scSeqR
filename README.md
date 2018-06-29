# scSeqR (Version: 0.99.0)

Maintainers: Alireza Khodadadi-Jamayran and Aristotelis Tsirigos.

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
     norm.method = "ranked.glsf",
     top.rank = 500) # best for scRNA-Seq

# more examples
#my.obj <- norm.data(my.obj, norm.method = "global.glsf") # best for bulk RNA-Seq 
#my.obj <- norm.data(my.obj, norm.method = "rpm", rpm.factor = 100000) # best for bulk RNA-Seq
#my.obj <- norm.data(my.obj, norm.method = "spike.in", spike.in.factors = NULL)
#my.obj <- norm.data(my.obj, norm.method = "no.norm") # if the data is already normalized
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

scSeqR let's you perform 2 dimentional and 3 dementional PCA and tSNE clustering (recomanded). You also have the option of distance based clustering using "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski" methods.

```r
# tSNE
my.obj <- cluster.data(my.obj, clust.method = "gene.model", gene.list = "my_model_genes.txt", clust.type = "tsne")

# PCA
my.obj <- cluster.data(my.obj, clust.method = "gene.model", gene.list = "my_model_genes.txt", clust.type = "pca")
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
	plot.type = "pca",
	col.by = "conditions",
	clust.dim = 2,
	interactive = F)
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/tSNE_conds.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/PCA_conds.png" width="400"/>      
</p>

- Find optimal number of clusters

scSeqR allows to choose from 3 different algorisms ("elbow.wss", "silhouette", "gap.stat") to find the optimal number of clusters. 

```r
opt.clust.num(my.obj, max.clust = 10, clust.type = "tsne", opt.method = "elbow.wss")
opt.clust.num(my.obj, max.clust = 10, clust.type = "tsne", opt.method = "silhouette")
# opt.clust.num(my.obj, max.clust = 10, gap.stat.nboot = 200, clust.type = "tsne",opt.method = "gap.stat")
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/optim_clust_num1.png" width="800" />
</p>
   
- Assign clusters

```r
# tSNE
my.obj <- assign.clust(my.obj, clust.num = 4,clust.type = "tsne")

# PCA
my.obj <- assign.clust(my.obj, clust.num = 4,clust.type = "pca")
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
	interactive = F,
	density = F)  
	
# more examples 
cluster.plot(my.obj,plot.type = "tsne",col.by = "clusters",clust.dim = 3,interactive = F)
cluster.plot(my.obj,plot.type = "pca",col.by = "clusters",clust.dim = 2,interactive = F) 
cluster.plot(my.obj,plot.type = "tsne",col.by = "conditions",clust.dim = 2,interactive = F,density = T)	
cluster.plot(my.obj,plot.type = "tsne",col.by = "clusters",clust.dim = 2,interactive = F,density = T)	
```

To see the above made interactive plots click on these links: [2Dplot](https://rawgit.com/rezakj/scSeqR/master/doc/tSNE_plot2d_clustered.html) and [3Dplot](https://rawgit.com/rezakj/scSeqR/dev/doc/tSNE_3D_clusters.html)
        
<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/tSNE_2D_clusters.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/PCA_2D_clusters.png" width="400"/> 
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/tSNE_conds_density.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/tSNE_density.png" width="400"/> 	
</p>

        
- Avrage expression per cluster

```r
my.obj <- clust.avg.exp(my.obj, clust.type = "tsne", clust.dim = 2)
head(my.obj@clust.avg)
#      gene  cluster_1  cluster_2    cluster_3   cluster_4
#1     A1BG 0.07139890 0.04338571 0.0676852530 0.076330913
#2 A1BG.AS1 0.01291263 0.01490700 0.0107435934 0.003482156
#3     A1CF 0.00000000 0.00000000 0.0000000000 0.000000000
#4      A2M 0.00247426 0.00000000 0.0009037397 0.001937402
#5  A2M.AS1 0.01037358 0.00000000 0.0044210665 0.041431989
#6    A2ML1 0.00000000 0.00000000 0.0000000000 0.000000000
```

- Save your object

```r
save(my.obj, file = "my.obj.Robj")
```        

- Find marker genes

```r
marker.genes <- find.markers(my.obj,
	data.type = "tsne",
	fold.change = 2,
	padjval = 0.1)

head(marker.genes)
#            baseMean    baseSD AvExpInCluster AvExpInOtherClusters foldChange
#LINC00176 0.06768702 0.2952005     0.15649818          0.010201374  15.340892
#ADTRP     0.04528819 0.2545610     0.10278008          0.008074863  12.728400
#TSHZ2     0.04839484 0.2631226     0.10975496          0.008677692  12.647944
#WNT7A     0.01022659 0.1024340     0.02252658          0.002265051   9.945287
#EPHX2     0.05590041 0.2518293     0.12250585          0.012788071   9.579697
#MAL       0.19356697 1.0878996     0.42139913          0.046095870   9.141798
#          log2FoldChange         pval         padj clusters      gene  cluster_1
#LINC00176       3.939311 8.339096e-25 1.870459e-21        1 LINC00176 0.15649818
#ADTRP           3.669979 8.418255e-15 1.867169e-11        1     ADTRP 0.10278008
#TSHZ2           3.660831 1.220141e-15 2.712374e-12        1     TSHZ2 0.10975496
#WNT7A           3.314013 3.441443e-05 7.344040e-02        1     WNT7A 0.02252658
#EPHX2           3.259980 2.914413e-20 6.522456e-17        1     EPHX2 0.12250585
#MAL             3.192478 1.672668e-12 3.698268e-09        1       MAL 0.42139913
#            cluster_2   cluster_3   cluster_4
#LINC00176 0.003401901 0.003570966 0.022241469
#ADTRP     0.007460135 0.010882054 0.005106448
#TSHZ2     0.002285751 0.008833551 0.012387695
#WNT7A     0.000000000 0.000000000 0.006342935
#EPHX2     0.000000000 0.005966443 0.028705683
#MAL       0.026124748 0.018165270 0.091529773
```

- Plot genes

```r
# tSNE 2D
gene.plot(my.obj, gene = "MS4A1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "Cebpb_scatter_plot")
# PCA 2D	
gene.plot(my.obj, gene = "MS4A1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "Cebpb_scatter_plot",
	plot.data.type = "pca")
	
# tSNE 3D	
gene.plot(my.obj, gene = "MS4A1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "Cebpb_scatter_plot",
	clust.dim = 3)
	
# Box Plot
gene.plot(my.obj, gene = "MS4A1", 
	box.to.test = 0, 
	box.pval = "sig.values",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")
	
# Bar plot (to visualize fold changes)	
gene.plot(my.obj, gene = "MS4A1", 
	col.by = "clusters",
	plot.type = "barplot",
	interactive = F,
	out.name = "bar_plot")
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/MS4A1_tSNE.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/MS4A1_PCA.png" width="400"/> 
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/MS4A1_box.png" width="400"/> 
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/MS4A1_bar.png" width="400"/>
</p>

- Heatmap

```r
MyGenes <- top.markers(marker.genes, topde = 20, SDmin = 5)

heatmap.plot (my.obj,
	gene = MyGenes,
	plot.data.type = "tsne",
	clust.dim = 2,
	cluster.by = "clusters",
	cluster.rows = T,
	heat.colors = c("yellow" ,"black", "green"))
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/heatmap.png" width="800" height="800" />
</p>


- Differential Expression Analysis 

```r
diff.res <- diff.exp(my.obj, de.by = "clusters", cond.1 = c(1,4), cond.2 = c(2))
head(diff.res)
#            baseMean         1_4          2 foldChange log2FoldChange         pval
#A1BG     0.067851093 0.073152434 0.04338571  0.5930863     -0.7536861 4.899112e-02
#A1BG.AS1 0.010512038 0.009559705 0.01490700  1.5593579      0.6409521 4.978422e-01
#A1CF     0.000000000 0.000000000 0.00000000        NaN            NaN          NaN
#A2M      0.001876722 0.002283384 0.00000000  0.0000000           -Inf 4.845531e-02
#A2M.AS1  0.017602004 0.021416137 0.00000000  0.0000000           -Inf 3.324006e-06
#A2ML1    0.000000000 0.000000000 0.00000000        NaN            NaN          NaN
#               padj
#A1BG     1.00000000
#A1BG.AS1 1.00000000
#A1CF            NaN
#A2M      1.00000000
#A2M.AS1  0.04989001
#A2ML1           NaN

# more examples 
diff.res <- diff.exp(my.obj, de.by = "conditions", cond.1 = c("KO"), cond.2 = c("WT"))
```

- Volcano and MA plots 

```r
# Volcano Plot 
volcano.ma.plot(diff.res,
	sig.value = "pval",
	sig.line = 0.05,
	plot.type = "volcano",
	interactive = F)

# MA Plot
volcano.ma.plot(diff.res,
	sig.value = "pval",
	sig.line = 0.05,
	plot.type = "ma",
	interactive = F)
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/volc_plot.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/MA_plot.png" width="400"/>      
</p>

