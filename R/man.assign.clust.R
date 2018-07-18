#' Assign the optimal number of clusters to hierarchical clustering.
#'
#' This function takes an object of class scSeqR and finds optimal number of clusters based on three methods.
#' @param x An object of class scSeqR.
#' @param clust.num Number of clusters
#' @param  clust.type Choose from "tsne","pca" or "distance", defult = "tsne".
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' assign.clust(my.obj, clust.num = 7)
#' }
#' @export
man.assign.clust <- function (x = NULL,
                          clust.num = 0) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  if ( clust.num == 0) {
    stop("please provite the optimal number of clusters")
  }
  ClustNum = clust.num
#####
    DATA <- (x@tsne.data)
    fit_cluster_hierarchical = hclust(dist(scale(DATA)))
    DATA$clusters = factor(cutree(fit_cluster_hierarchical, k=ClustNum))
    DATAclusters = as.data.frame(cutree(fit_cluster_hierarchical, k = ClustNum))
    colnames(DATAclusters) <- c("clusters")
    attributes(x)$best.clust <- DATAclusters
  #
  return(x)
}
