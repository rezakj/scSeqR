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
assign.clust <- function (x = NULL,
                          clust.num = 0,
                          clust.type = "tsne") {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  if ( clust.num == 0) {
    stop("please provite the optimal number of clusters")
  }
  ClustNum = clust.num
  #  if (clust.dim == 3) {
  if (clust.type == "tsne") {
    DATA <- (x@tsne.data)[1:2]
    fit_cluster_hierarchical = hclust(dist(scale(DATA)))
    DATA$clusters = factor(cutree(fit_cluster_hierarchical, k = ClustNum))
    DATAclusters = as.data.frame(cutree(fit_cluster_hierarchical, k = ClustNum))
    colnames(DATAclusters) <- c("clusters")
    tsne.data <- DATA
    attributes(x)$tsne.data <- tsne.data
    # merge with 3 d
    DATA <- (x@tsne.data.3d)[1:3]
    mymrgd <- cbind(DATA, tsne.data[3])
    attributes(x)$tsne.data.3d <- mymrgd
  }
  if (clust.type == "pca") {
    DATA <- (x@pca.data)[1:2]
    fit_cluster_hierarchical = hclust(dist(scale(DATA)))
    DATA$clusters = factor(cutree(fit_cluster_hierarchical, k=ClustNum))
    DATAclusters = as.data.frame(cutree(fit_cluster_hierarchical, k = ClustNum))
    colnames(DATAclusters) <- c("clusters")
    PCAdata <- DATA
    attributes(x)$pca.data <- PCAdata
    # merge with 3 d
    DATA <- (x@pca.data.3d)[1:3]
    mymrgd <- cbind(DATA, PCAdata[3])
    attributes(x)$pca.data.3d <- mymrgd
  }
  #  }
  if (clust.type == "distance") {
    DATA <- x@dist.data
    fit_cluster_hierarchical = hclust(dist(scale(DATA)))
    DATA$clusters = factor(cutree(fit_cluster_hierarchical, k=ClustNum))
    dist.data <- DATA
    attributes(x)$dist.data <- dist.data
  }
  return(x)
}
