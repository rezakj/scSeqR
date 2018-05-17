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
#' filter.cells(my.obj,
#'    min.mito = 0,
#'    max.mito = 1,
#'    min.genes = 0,
#'    max.genes = Inf,
#'    min.umis = 0,
#'    max.umis = Inf)
#' }
#' @export
filter.cells <- function (x = NULL, min.mito = 0, max.mito = 1, min.genes = 0, max.genes = Inf, min.umis = 0, max.umis = Inf) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  DATA <- x@raw.data
  MAXmito <- as.character(subset(x@stats, x@stats$mito.percent > max.mito)$CellIds)
  MINmito <- as.character(subset(x@stats, x@stats$mito.percent < min.mito)$CellIds)
  MAXbad <- as.character(subset(x@stats, x@stats$nGenes > max.genes)$CellIds)
  MINbad <- as.character(subset(x@stats, x@stats$nGenes < min.genes)$CellIds)
  MAXbadUMI <- as.character(subset(x@stats, x@stats$UMIs > max.umis)$CellIds)
  MINbadUMI <- as.character(subset(x@stats, x@stats$UMIs < min.umis)$CellIds)
  BadMito <- c(MAXmito,MINmito)
  BadGenes <- c(MAXbad,MINbad)
  BadUMIs <- c(MAXbadUMI,MINbadUMI)
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
  attributes(x)$main.data <- DATA
  return(x)
}
