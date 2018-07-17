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
reassign.clust <- function (x = NULL,
                          change.clust = 0,
                          to.clust = 0,
                          clust.type = "tsne") {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
###################
  if (clust.type == "tsne") {
    DATA <- (x@tsne.data)
    DATA$clusters <- gsub(change.clust,to.clust,DATA$clusters)
    attributes(x)$tsne.data <- DATA
    # 3d
    DATA <- (x@tsne.data.3d)
    DATA$clusters <- gsub(change.clust,to.clust,DATA$clusters)
    attributes(x)$tsne.data.3d <- DATA
  }
  if (clust.type == "pca") {
    DATA <- (x@pca.data)
    DATA$clusters <- gsub(change.clust,to.clust,DATA$clusters)
    attributes(x)$pca.data <- DATA
    # merge with 2 d
    DATA <- (x@pca.data.3d)
    DATA$clusters <- gsub(change.clust,to.clust,DATA$clusters)
    attributes(x)$pca.data.3d <- DATA
  }
##############
  return(x)
}
