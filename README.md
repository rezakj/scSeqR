# scSeqR v0.99.0

Authors: Alireza Khodadadi-Jamayran and Aristotelis Tsirigos.

### Single Cell Sequencing R package (scSeqR)

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/out.gif" width="700" height="600" />
</p>

scSeqR is an R package that can analyze single cell sequencing data types (i.e [scRNA-seq](https://en.wikipedia.org/wiki/Single_cell_sequencing#Single-cell_RNA_sequencing_(scRNA-seq))) and large numeric [matrix](https://en.wikipedia.org/wiki/Matrix_(mathematics)) 
files (i.e. count tables with many samples from [TCGA](https://cancergenome.nih.gov/)). The program inputs single cell data in [10X format](https://www.10xgenomics.com/), large numeric **matrix files** and **data frames** and helps you to perform QC, filtering, visualization, normalization, clustering, differential expression analysis and find positive and negative markers for each cluster. scSeqR, allows you to choose from **multiple normalization** methods and **spike-in normalization** depending on your data type. Alternatively, you can also use 
**already normalized** data. You also have the option of choosing from a variaty of clustering algorithms. 

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
     "gridExtra",
     "scatterplot3d",
     "RColorBrewer",
     "pheatmap"))
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

## To see interactive version of the plots click on these links: [mito-UMIs plot](https://rawgit.com/rezakj/scSeqR/dev/doc/mito-umi-plot.html) and [gene-UMIs plot](https://rawgit.com/rezakj/scSeqR/dev/doc/gene-umi-plot.html)

- Filter cells. 

scSeqR allows you to filter based on library sizes (UMIs), number of genes per cell, percent mitochondial content and even based on one or more genes, or cell ids. 

Befor filtering let's do a quick test and see how the data looks like. This might be very helpful to check the gene or genes of your interest and see in howmany cells they are expressed in. This could also be helful to choose the genes for cell sorting or filteering that are needed in some studies. 

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

Let's say that you are only interested in the cells in which RPL13 "OR" RPL10 are expressed, while have a maximum mito percent of 0.05 and a minimum number of 250 genes expressed and a maximum of 2400 genes. Maybe you also want to filter a cell with this id: WT_AAACATACAACCAC.1 

```r
my.obj <- cell.filter(my.obj,
	filter.by.gene = c("RPL13","RPL10"),
	filter.by.cell.id = c("WT_AAACATACAACCAC.1"),
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


- Perform PCA and tSNE

```r
# PCA
my.obj <- run.pca(my.obj, clust.method = "gene.model", gene.list = "my_model_genes.txt")

# tSNE
my.obj <- run.tsne(my.obj, clust.method = "gene.model", gene.list = "my_model_genes.txt")
```        

- Cluster the data

Here we cluster the first 10 dimensions of the data which is converted to principal components, to do this, you have the option of clustering your data based on the following methods: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans"

 For the distance calculation used for clustering, you have the following options: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski" or "NULL"

 With the following indexing methods: "kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw"

We recomand to use the defult options as below:

```r
my.obj <- run.clustering(my.obj, 
	clust.method = "kmeans", 
	dist.method = "euclidean",
	index.method = "silhouette",
	max.clust = 20,
	dims = 1:10)

# number of clusters found and assigned

my.obj@cluster.data$Best.nc
#Number_clusters     Value_Index 
#         7.0000          0.2849 
```

 - Optional manual clustering or renaming the clusters 
 
 You also have the option of manual hirarchical clustering or renaming the clusters. It is highly recomanded to not use this method as the above method is much more accurate. 
To do this you might need to see what is the optimal number of clusters. 

```r
##### Find optimal number of clusters for hierarchical clustering
#opt.clust.num(my.obj, max.clust = 10, clust.type = "tsne", opt.method = "silhouette")
##### Manual clustering 
#my.obj <- man.assign.clust(my.obj, clust.num = 7)
##### re-assign clusters 
#my.obj <- change.clust(my.obj,change.clust = 1,to.clust = 20)
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/optim_clust_num1.png" width="800" />
</p>

- Visualize conditions

As we artificially made 3 conditions by randomly dividing the sample into 3. All the conditions should be looking similar.

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

- Visualize clusters

```r
# 2D
cluster.plot(my.obj,
	cell.size = 1,
	plot.type = "tsne",
	cell.color = "black",
	back.col = "white",
	col.by = "clusters",
	cell.transparency = 0.5,
	clust.dim = 2,
	interactive = F)

# 3D
cluster.plot(my.obj,
	plot.type = "tsne",
	col.by = "clusters",
	clust.dim = 3,
	interactive = F,
	density = F,
	angle = 100)
	
# intractive 2D
cluster.plot(my.obj,
	plot.type = "tsne",
	col.by = "clusters",
	clust.dim = 2,
	interactive = T,
	out.name = "tSNE_2D_clusters")

# intractive 3D
cluster.plot(my.obj,
	plot.type = "tsne",
	col.by = "clusters",
	clust.dim = 3,
	interactive = T,
	out.name = "tSNE_3D_clusters")

# Density plot for clusters 
cluster.plot(my.obj,
	plot.type = "pca",
	col.by = "clusters",
	interactive = F,
	density=T)

# Density plot for conditions 
cluster.plot(my.obj,
	plot.type = "pca",
	col.by = "conditions",
	interactive = F,
	density=T)
```
## To see the above made interactive plots click on these links: [2Dplot](https://rawgit.com/rezakj/scSeqR/dev/doc/tSNE_2D_clusters.html) and [3Dplot](https://rawgit.com/rezakj/scSeqR/dev/doc/tSNE_3D_clusters.html)
        
<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/tSNE_2D_clusters.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/tSNE_3D.png" width="400"/> 
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/density_conditions.png" width="400"/>
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/density_clusters.png" width="400"/> 	
</p>

        
- Avrage expression per cluster

```r
my.obj <- clust.avg.exp(my.obj)

head(my.obj@clust.avg)
#     gene   cluster_1   cluster_2   cluster_3 cluster_4   cluster_5   cluster_6  cluster_7
#1     A1BG 0.074805398 0.083831677 0.027234682         0 0.088718322 0.026671084 0.04459271
#2 A1BG.AS1 0.013082859 0.012882983 0.005705715         0 0.003077574 0.000000000 0.01498637
#3     A1CF 0.000000000 0.000000000 0.000000000         0 0.000000000 0.000000000 0.00000000
#4      A2M 0.002350504 0.000000000 0.003284837         0 0.000000000 0.006868043 0.00000000
#5  A2M.AS1 0.009734684 0.006208601 0.000000000         0 0.041558965 0.055534823 0.00000000
#6    A2ML1 0.000000000 0.000000000 0.000000000         0 0.000000000 0.000000000 0.00000000
```

- Save your object

```r
save(my.obj, file = "my.obj.Robj")
```        

- Find marker genes

```r
marker.genes <- find.markers(my.obj,
	fold.change = 2,
	padjval = 0.1)

dim(marker.genes)
# [1] 1070   17

head(marker.genes)
#             baseMean     baseSD AvExpInCluster AvExpInOtherClusters foldChange log2FoldChange
#WNT7A     0.010229718 0.10242607     0.02380399         0.0006494299   36.65368       5.195886
#KRT1      0.020771859 0.18733189     0.04765674         0.0017973688   26.51473       4.728722
#TSHZ2     0.048409666 0.26309988     0.11075764         0.0044064635   25.13527       4.651641
#ANKRD55   0.008418143 0.09648853     0.01920420         0.0008056872   23.83580       4.575058
#LINC00176 0.067707756 0.29517298     0.15445793         0.0064822618   23.82778       4.574573
#MAL       0.193626263 1.08780664     0.42990647         0.0268672596   16.00113       4.000102
#                  pval         padj clusters      gene  cluster_1   cluster_2   cluster_3 cluster_4
#WNT7A     1.258875e-06 2.772043e-03        1     WNT7A 0.02380399 0.000000000 0.000000000         0
#KRT1      1.612082e-07 3.577209e-04        1      KRT1 0.04765674 0.000000000 0.000000000         0
#TSHZ2     3.647481e-18 8.341790e-15        1     TSHZ2 0.11075764 0.009321206 0.007981973         0
#ANKRD55   3.871894e-05 8.409755e-02        1   ANKRD55 0.01920420 0.000000000 0.000000000         0
#LINC00176 2.439482e-27 5.620565e-24        1 LINC00176 0.15445793 0.003689827 0.003429306         0
#MAL       2.417583e-15 5.500001e-12        1       MAL 0.42990647 0.021592215 0.010139857         0
#            cluster_5  cluster_6   cluster_7
#WNT7A     0.002806920 0.00000000 0.000000000
#KRT1      0.005251759 0.00000000 0.002596711
#TSHZ2     0.000000000 0.00000000 0.002297921
#ANKRD55   0.003482284 0.00000000 0.000000000
#LINC00176 0.013798598 0.00910279 0.003420015
#MAL       0.041585770 0.03214897 0.026263849
```

- Plot genes

```r
# tSNE 2D
gene.plot(my.obj, gene = "MS4A1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot")
# PCA 2D	
gene.plot(my.obj, gene = "MS4A1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
	plot.data.type = "pca")
	
# tSNE 3D	
gene.plot(my.obj, gene = "MS4A1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "scatter_plot",
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


- Multiple plots

```r
genelist = c("PPBP","LYZ","MS4A1","GNLY","LTB","NKG7","IFITM2","CD14","S100A9")
###
library(gridExtra)
for(i in genelist){
	MyPlot <- gene.plot(my.obj, gene = i, 
		plot.type = "scatterplot",
		interactive = F,
		out.name = "Cebpb_scatter_plot")
	eval(call("<-", as.name(i), MyPlot))
}
### plot 
grid.arrange(PPBP,LYZ,MS4A1,GNLY,LTB,NKG7,IFITM2,CD14,S100A9)
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/list1.png" width="800" height="800" />
</p>


- Heatmap

```r
# find top genes
MyGenes <- top.markers(marker.genes, topde = 10, min.base.mean = 0.8)
# plot
heatmap.plot (my.obj, gene = MyGenes)
```

<p align="center">
  <img src="https://github.com/rezakj/scSeqR/blob/dev/doc/heatmap.png" width="800" height="800" />
</p>


- Differential Expression Analysis 

```r
diff.res <- diff.exp(my.obj, de.by = "conditions", cond.1 = c("WT"), cond.2 = c("KO"))
head(diff.res)
#            baseMean          WT          KO foldChange log2FoldChange       pval padj
#A1BG     0.068250244 0.073606778 0.062875387  0.8542065     -0.2273433 0.46284942    1
#A1BG.AS1 0.009902209 0.012277809 0.007518483  0.6123636     -0.7075395 0.31076602    1
#A1CF     0.000000000 0.000000000 0.000000000        NaN            NaN        NaN  NaN
#A2M      0.001555311 0.003105320 0.000000000  0.0000000           -Inf 0.09170874    1
#A2M.AS1  0.010411651 0.005819828 0.015019181  2.5806916      1.3677577 0.13272473    1
#A2ML1    0.000000000 0.000000000 0.000000000        NaN            NaN        NaN  NaN

diff.res <- diff.exp(my.obj, de.by = "clusters", cond.1 = c(1,4), cond.2 = c(2))
head(diff.res)
#            baseMean         1_4          2 foldChange log2FoldChange       pval padj
#A1BG     0.073234544 0.094378843 0.02684924  0.2844837     -1.8135823 0.01472597    1
#A1BG.AS1 0.002151004 0.003131519 0.00000000  0.0000000           -Inf 0.31800136    1
#A1CF     0.000000000 0.000000000 0.00000000        NaN            NaN        NaN  NaN
#A2M      0.002164828 0.000000000 0.00691392        Inf            Inf 0.31882994    1
#A2M.AS1  0.042598420 0.036532388 0.05590578  1.5303072      0.6138213 0.52322284    1
#A2ML1    0.000000000 0.000000000 0.00000000        NaN            NaN        NaN  NaN

diff.res <- diff.exp(my.obj, de.by = "clustBase.condComp", cond.1 = c("WT"), cond.2 = c("KO"), base.cond = 1)
head(diff.res)
#            baseMean WT.inClust.1 KO.inClust.1 foldChange log2FoldChange      pval padj
#A1BG     0.085041487  0.041275761   0.12841293   3.111098       1.637424 0.1363094    1
#A1BG.AS1 0.004973589  0.009992392   0.00000000   0.000000           -Inf 0.3195253    1
#A1CF     0.000000000  0.000000000   0.00000000        NaN            NaN       NaN  NaN
#A2M      0.000000000  0.000000000   0.00000000        NaN            NaN       NaN  NaN
#A2M.AS1  0.028292066  0.011177204   0.04525274   4.048664       2.017446 0.2468776    1
#A2ML1    0.000000000  0.000000000   0.00000000        NaN            NaN       NaN  NaN

diff.res <- diff.exp(my.obj, de.by = "condBase.clustComp", cond.1 = c(1), cond.2 = c(2), base.cond = "WT")
head(diff.res)
#            baseMean 1.inCond.WT 2.inCond.WT foldChange log2FoldChange      pval padj
#A1BG     0.042345533 0.041275761  0.04452470   1.078713      0.1093111 0.9314472    1
#A1BG.AS1 0.006702214 0.009992392  0.00000000   0.000000           -Inf 0.3195253    1
#A1CF     0.000000000 0.000000000  0.00000000        NaN            NaN       NaN  NaN
#A2M      0.006745287 0.000000000  0.02048569        Inf            Inf 0.3218544    1
#A2M.AS1  0.007496905 0.011177204  0.00000000   0.000000           -Inf 0.3195253    1
#A2ML1    0.000000000 0.000000000  0.00000000        NaN            NaN       NaN  NaN
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

