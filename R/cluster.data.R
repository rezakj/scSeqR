#' Cluster based on tSNE
#'
#' This function takes an object of class scSeqR and clusters the normalized data.
#' @param x An object of class scSeqR.
#' @param clust.method a method to choose your genes for clustering. There are three options "base.mean.rank","dispersed.genes","both" defults is "base.mean.rank".
#' @param top.rank number of ranked genes you want to use for clustering, defult = 500.
#' @param clust.dim number of dimentions you want to consider for clustering, defult =2.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' cluster.data(my.obj, clust.method = "base.mean.rank", top.rank = 500, clust.dim = 3)
#' }
#' @import Rtsne
#' @export
cluster.data <- function (x = NULL, clust.method = "base.mean.rank", top.rank = 500, clust.dim = 2) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  if (clust.method == "base.mean.rank") {
    DATA <- x@main.data
    raw.data.order <- DATA[ order(rowMeans(DATA), decreasing = T), ]
    topGenes = head(raw.data.order,top.rank)
    NormLog <- log(topGenes+0.1)
    TopNormLogScale <- as.data.frame(scale(NormLog))
    TransPosed <- t(TopNormLogScale)
    tsne <- Rtsne(TransPosed, check_duplicates = FALSE, dims = clust.dim)
    tsne.data = as.data.frame(tsne$Y)
    tsne.data = cbind(cells = row.names(TransPosed),tsne.data)
  }
  attributes(x)$tsne.data <- tsne.data
  return(x)
}

