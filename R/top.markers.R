#' Load 10X data as data.frame
#'
#' This function takes 10X data files barcodes.tsv, genes.tsv and matrix.mtx and converts them to proper matrix file for scSeqR.
#' @param dir.10x A directory that includes the 10X barcodes.tsv, genes.tsv and matrix.mtx files.
#' @param gene.name Should be either geneSymbol or ensembleID.
#' @return The data frame object
#' @examples
#' \dontrun{
#' load10x("/hg19", gene.name = "geneSymbol")
#' }
#' @import Matrix
#' @export
top.markers <- function (x = NULL, topde = 5, SDmin = 1) {
  MyClusts <- (unique(x$clusters))
  for (i in MyClusts) {
    DATA <- subset(x, x$clusters == i)
    DATA <- subset(DATA, DATA$baseSD > SDmin)
    DATA <- as.character(head(DATA,topde)$gene)
    DatNmaes=paste("topgenes",i,sep="_")
    eval(call("<-", as.name(DatNmaes), DATA))
  }
  # cat them
  filenames <- ls(pattern="topgenes_")
  datalist <- mget(filenames)
  topGenes <- as.character(do.call("c", datalist))
    return(topGenes)
  }


