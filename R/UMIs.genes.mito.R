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
  UMIs <- colSums(x)
  UMIs.df <- as.data.frame(UMIs)
  nGenes <- sapply(x, function(x) length(as.numeric(subset(x, x != 0))))
  nGenes.df <- as.data.frame(nGenes)
  mito.genes <- grep(pattern = "^mt\\-", x = rownames(x), value = TRUE, ignore.case = TRUE)
  mito <- subset(x,rownames(x) %in% mito.genes)
  mitoSiz <- colSums(mito)
  mito.percent <- (mitoSiz/libSiz)
  mito.percent.df <- as.data.frame(mito.percent)
  return(c(mito.percent.df,nGenes.df,UMIs.df))
}
