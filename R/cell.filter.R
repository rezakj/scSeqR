#' Filter cells
#'
#' This function takes an object of class scSeqR and filters the raw data based on the number of UMIs, genes per cell and percentage of mitochondrial genes per cell.
#' @param x An object of class scSeqR.
#' @param min.mito Min ratio of mit genes, defult = 0.
#' @param max.mito Max ratio of mit genes, defult = 1.
#' @param min.genes Min number genes per cell, defult = 0.
#' @param max.genes Max number genes per cell, defult = Inf.
#' @param min.umis Min number UMIs per cell, defult = 0.
#' @param max.umis Max number UMIs per cell, defult = Inf.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' cell.filter(my.obj,
#'    min.mito = 0,
#'    max.mito = 1,
#'    min.genes = 0,
#'    max.genes = Inf,
#'    min.umis = 0,
#'    max.umis = Inf)
#' }
#'
#' @export
cell.filter <- function (x = NULL,
                         min.mito = 0,
                         max.mito = 1,
                         min.genes = 0,
                         max.genes = Inf,
                         min.umis = 0,
                         max.umis = Inf,
                         filter.by.gene = "character",
                         filter.by.gene.min = 1) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
# take input from
  DATA <- x@raw.data
# filter by gene
  if (filter.by.gene[1] != "character"){
    genExist <- dim(subset(DATA, row.names(DATA) %in% filter.by.gene))[1]
    numberOfgenes <- length(filter.by.gene)
    if (genExist != numberOfgenes) {
      stop("your data lacks the gene/genes in filter.by.gene.")
    }
  gene.filt <- t(subset(DATA, row.names(DATA) %in% filter.by.gene) >= filter.by.gene.min)
#  gene.filted <- as.data.frame(rowSums(gene.filt) > 0L)
  gene.filted <- as.data.frame(apply(gene.filt, 1, any))
  colnames(gene.filted) <- "gene.filter"
  gene.filt.cells <- row.names(subset(gene.filted, gene.filter == F))
  goneCells <- length(gene.filt.cells)
  allCells <- length(colnames(DATA))
  }
# filter by mito
  MAXmito <- as.character(subset(x@stats, x@stats$mito.percent > max.mito)$CellIds)
  MINmito <- as.character(subset(x@stats, x@stats$mito.percent < min.mito)$CellIds)
# filter by nGene
  MAXbad <- as.character(subset(x@stats, x@stats$nGenes > max.genes)$CellIds)
  MINbad <- as.character(subset(x@stats, x@stats$nGenes < min.genes)$CellIds)
# filter by UMI
  MAXbadUMI <- as.character(subset(x@stats, x@stats$UMIs > max.umis)$CellIds)
  MINbadUMI <- as.character(subset(x@stats, x@stats$UMIs < min.umis)$CellIds)
# take all bad genes
  BadMito <- c(MAXmito,MINmito)
  BadGenes <- c(MAXbad,MINbad)
  BadUMIs <- c(MAXbadUMI,MINbadUMI)
# filter
  if (length(BadMito) == 0) {
    print("No mito filter")
  } else {
    DATA <- DATA[ , -which(names(DATA) %in% BadMito)]
    print(paste("cells with min mito ratio of",min.mito,"and max mito ratio of",max.mito,"were filtered."))
  }
  if (length(BadGenes) == 0) {
    print("No gene number filter")
  } else {
    DATA <- DATA[ , -which(names(DATA) %in% BadGenes)]
    print(paste("cells with min genes of",min.genes,"and max genes of",max.genes,"were filtered."))
  }
  if (length(BadUMIs) == 0) {
    print("No UMI number filter")
  } else {
    DATA <- DATA[ , -which(names(DATA) %in% BadUMIs)]
    print(paste("cells with min UMIs of",min.umis,"and max UMIs of",max.umis,"were filtered."))
  }
  if (!exists("gene.filt.cells")) {
    print("No cell filter by provided gene/genes")
  } else {
    DATA <- DATA[ , -which(names(DATA) %in% gene.filt.cells)]
    GenesForFilter = paste(filter.by.gene, collapse=",")
    print(paste(goneCells,"cells out of original",allCells,
                " cells were filtered out as their expression was less than",
                filter.by.gene.min,"for",GenesForFilter))
  }
# make a filter parameters file
  GenesForFilter = paste(filter.by.gene, collapse=",")
  FilterFile <- paste("Filters are:
  - min mito",min.mito,"max mito",max.mito,"
  - min # of genes",min.genes,"max # of genes",max.genes,"
  - min UMIs",min.umis,"max UMIs",max.umis,"
  - genes filtered by
  ",GenesForFilter)
# write filter parameters
  write.table((FilterFile),file="filters_set.txt", row.names =F, quote = F, col.names = F)
  print("filters_set.txt file has beed generated and includes the filters set for this experiment.")
# return data
  attributes(x)$main.data <- DATA
  return(x)
}
