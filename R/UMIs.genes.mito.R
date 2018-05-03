#' Calculate the number of UMIs and genes per cell and percentage of mitochondrial genes per cell.
#'
#' This function takes data frame and calculates the number of UMIs, genes per cell and percentage of mitochondrial genes per cell.
#' @param x A data frame containing gene counts for cells.
#' @return The data frame object
#' @examples
#' \dontrun{
#' UMIs.genes.mit(my.data)
#' }
#' @export
UMIs.genes.mito <- function (x = NULL) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  DATA <- x@raw.data
  UMIs <- colSums(DATA)
  nGenes <- sapply(DATA, function(DATA) length(as.numeric(subset(DATA, DATA != 0))))
  mito.genes <- grep(pattern = "^mt\\-", x = rownames(DATA), value = TRUE, ignore.case = TRUE)
  mito <- subset(DATA,rownames(DATA) %in% mito.genes)
  mitoSiz <- colSums(mito)
  mito.percent <- (mitoSiz/UMIs)
  QC <- list(colnames(DATA),as.numeric(nGenes), as.numeric(UMIs),as.numeric(mito.percent))
  names(QC) <- c("CellIds","nGenes","UMIs","mito.percent")
  STATS <- as.data.frame(QC)
  attributes(x)$stats <- STATS
  return(x)
}
