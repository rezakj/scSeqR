#' Change the cluster number or re-name them
#'
#' This function re-names the clusters in the best.clust slot of the scSeqR object.
#' @param x An object of class scSeqR.
#' @param change.clust The name of the cluster to be changed.
#' @param  to.clust The new name for the cluster.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' my.obj <- change.clust(my.obj, change.clust = 3, to.clust = 1)
#' my.obj <- change.clust(my.obj, change.clust = 2, to.clust = "B Cell")
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
