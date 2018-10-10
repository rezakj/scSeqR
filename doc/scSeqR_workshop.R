# https://genome.med.nyu.edu/results/external/temp/scSeqR_test/
# How to install scSeqR
install.packages("devtools")
library(devtools)
install_github("rezakj/scSeqR")


# set your working directory 
setwd("/your/download/directory")

####### download sample 10X data 
# save the URL as an object
sample.file.url = "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
# download the file
download.file(url = sample.file.url, 
     destfile = "pbmc3k_filtered_gene_bc_matrices.tar.gz", 
     method = "auto")  
# unzip the file. 
untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")

############################################# How to use scSeqR
library("scSeqR")
my.data <- load10x("filtered_gene_bc_matrices/hg19/",gene.name = "geneSymbol")

## Divide your sample into 3 as if you have 3 samples 
## This is to test scSeqR as if you have 3 conditions 
my.data <- load10x("../filtered_gene_bc_matrices/hg19/",gene.name = "geneSymbol")
head(my.data)[1:5]
dim(my.data)
  sample1 <- my.data[1:900]
  sample2 <- my.data[901:1800]
  sample3 <- my.data[1801:2700]
my.data <- data.aggregation(samples = c("sample1","sample2","sample3"), condition.names = c("WT","KO","Ctrl"))
#head(my.data)[1:5]

## Make a scSeqR object 
my.obj <- make.obj(my.data)

## see you object 
my.obj
## see the slots in the object 
slotNames(my.obj)

## QC
my.obj <- qc.stats(my.obj)

## visualize QC 
#png('myplot.png',width = 6, height = 6, units = 'in', res = 300)
stats.plot(my.obj,
	plot.type = "all.in.one",
	out.name = "UMI-plot",
	interactive = F,
	cell.color = "slategray3", 
	cell.size = 1, 
	cell.transparency = 0.5,
	box.color = "red",
	box.line.col = "green",
	back.col = "white")
# dev.off()
stats.plot(my.obj, plot.type = "point.mito.umi", interactive = F, out.name = "mito.umis.plot")

stats.plot(my.obj, plot.type = "point.gene.umi", interactive = F)

## filter the data 

my.obj <- cell.filter(my.obj,
	min.mito = 0,
	max.mito = 0.05,
	min.genes = 200,
	max.genes = 2400,
	min.umis = 0,
	max.umis = Inf)

## normalize your data for UMI differences 
my.obj <- norm.data(my.obj, norm.method = "ranked.glsf", top.rank = 500)

## scale your data 
my.obj <- data.scale(my.obj)

## repeat QC
my.obj <- qc.stats(my.obj,which.data = "main.data")

stats.plot(my.obj,
	plot.type = "all.in.one",
	out.name = "UMI-plot",
	interactive = F,
	cell.color = "slategray3", 
	cell.size = 1, 
	cell.transparency = 0.5,
	box.color = "red",
	box.line.col = "green",
	back.col = "white")

## QC and stats on genes 
my.obj <- gene.stats(my.obj, which.data = "main.data")

### make gene model
make.gene.model(my.obj, 
	dispersion.limit = 1.5, 
	base.mean.rank = 500, 
	no.mito.model = T, 
	mark.mito = T, 
	interactive = F,
	out.name = "gene.model")

## run PCA 
my.obj <- run.pca(my.obj, clust.method = "gene.model", gene.list = "my_model_genes.txt")

## significant PCs
opt.pcs.plot(my.obj)

## cluster the data 
my.obj <- run.clustering(my.obj, 
	clust.method = "ward.D", 
	dist.method = "euclidean",
	index.method = "kl",
	max.clust = 25,
	min.clust = 2,
	dims = 1:my.obj@opt.pcs)
## 

my.obj <- run.pc.tsne(my.obj, dims = 1:my.obj@opt.pcs)

# my.obj <- run.tsne(my.obj, clust.method = "gene.model", gene.list = "my_model_genes.txt")

## plot 
cluster.plot(my.obj,
	cell.size = 1,
	plot.type = "tsne",
	cell.color = "black",
	back.col = "white",
	col.by = "clusters",
	cell.transparency = 0.5,
	clust.dim = 2,
	interactive = F,
	density = F,
	out.name = "tSNE_3D_clusters")

clust.cond.info(my.obj, plot.type = "pie")

clust.cond.info(my.obj, plot.type = "bar")

## avrage per cluster 

my.obj <- clust.avg.exp(my.obj)
head(my.obj@clust.avg)

## QC

clust.stats.plot(my.obj, plot.type = "box.mito", interactive = F)
clust.stats.plot(my.obj, plot.type = "box.gene", interactive = F)
clust.stats.plot(my.obj, plot.type = "box.umi", interactive = F)

### save object 
save(my.obj, file = "my.obj.Robj")
 
# find marher genes

marker.genes <- find.markers(my.obj,
	fold.change = 2,
	padjval = 0.1)

head(marker.genes)

########## plot heatmap 

heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters")

## plot genes
### plot genes 
gene.plot(my.obj, gene = "MS4A1", 
	plot.type = "scatterplot",
	interactive = F,
	out.name = "myplot")

gene.plot(my.obj, gene = "MS4A1", 
	box.to.test = 0, 
	box.pval = "sig.values",
	col.by = "clusters",
	plot.type = "boxplot",
	interactive = F,
	out.name = "box_plot")

## multi plot 
genelist = c("PPBP","LYZ","MS4A1","GNLY","LTB","NKG7","IFITM2","CD14","S100A9")
noquote(paste(genelist, collapse = ","))

library(gridExtra)
for(i in genelist){
	MyPlot <- gene.plot(my.obj, gene = i, 
		plot.type = "scatterplot",
		interactive = F,
		out.name = "Cebpb_scatter_plot")
	eval(call("<-", as.name(i), MyPlot))
}

grid.arrange(PPBP,LYZ,MS4A1,GNLY,LTB,NKG7,IFITM2,CD14,S100A9)
rm(list = genelist)


### cell types 
Cluster = 2
MyGenes <- top.markers(marker.genes, topde = 40, min.base.mean = 0.2, cluster = Cluster)
#Name <- paste("ImmGen_heatmap_RNA_Cluster_",Cluster,".png",sep="")
#png(Name,width = 10, height = 10, units = 'in', res = 300)
imm.gen(immgen.data = "rna", gene = MyGenes, plot.type = "heatmap")
#dev.off()

#Name <- paste("ImmGen_pointPlot_RNA_Cluster_",Cluster,".png",sep="")
#png(Name,width = 10, height = 10, units = 'in', res = 300)
imm.gen(immgen.data = "rna", gene = MyGenes, plot.type = "point.plot")
#dev.off()


#Name <- paste("ImmGen_heatmap_ULI-RNA_Cluster_",Cluster,".png",sep="")
#png(Name,width = 20, height = 10, units = 'in', res = 300)
imm.gen(immgen.data = "uli.rna", gene = MyGenes, plot.type = "heatmap")
#dev.off()

#Name <- paste("ImmGen_pointPlot_ULI-RNA_Cluster_",Cluster,".png",sep="")
#png(Name,width = 10, height = 10, units = 'in', res = 300)
imm.gen(immgen.data = "uli.rna", gene = MyGenes, plot.type = "point.plot", top.cell.types = 50)
#dev.off()

## DE
diff.res <- diff.exp(my.obj, de.by = "clusters", cond.1 = c(1,4), cond.2 = c(2))
#diff.res <- diff.exp(my.obj, de.by = "conditions", cond.1 = c("WT"), cond.2 = c("KO"))
#diff.res <- diff.exp(my.obj, de.by = "clustBase.condComp", cond.1 = c("WT"), cond.2 = c("KO"), base.cond = 1)
#diff.res <- diff.exp(my.obj, de.by = "condBase.clustComp", cond.1 = c(1), cond.2 = c(2), base.cond = "WT")

volcano.ma.plot(diff.res,
	sig.value = "pval",
	sig.line = 0.05,
	plot.type = "volcano",
	interactive = F)

volcano.ma.plot(diff.res,
	sig.value = "pval",
	sig.line = 0.05,
	plot.type = "ma",
	interactive = F)

## renaming, removing and merging clusters 
my.obj <- change.clust(my.obj, change.clust = 3, to.clust = 2)
my.obj <- change.clust(my.obj, clust.reset = T)
my.obj <- change.clust(my.obj, change.clust = 7, to.clust = "B Cell")
my.obj <- clust.rm(my.obj, clust.to.rm = 1)























