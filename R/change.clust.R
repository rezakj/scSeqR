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
change.clust <- function (x = NULL,
                          change.clust = 0,
                          to.clust = 0) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
###################
    DATA <- (x@best.clust)
#    DATA$clusters <- gsub(change.clust,
#                          to.clust,
#                          DATA$clusters,
#                          fixed = T)
    DATA[DATA == change.clust] <- to.clust
    attributes(x)$best.clust <- DATA
##############
  return(x)
}
