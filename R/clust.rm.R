#' Create an object of class scSeqR.
#'
#' This function takes data frame and makes an object of class scSeqR.
#' @param x A data frame containing gene counts for cells.
#' @return An object of class scSeqR
#' @examples
#' \dontrun{
#' make.obj(my.data)
#' }
#'
#' @export
clust.rm <- function (x = NULL, clust.to.rm = "numeric") {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  if (clust.to.rm == "numeric") {
    stop("you should choose a cluster number to remove")
  }
  # find cluster ids to remove
  DATA <- x@best.clust
  clustersToGo <- row.names(subset(DATA, DATA$clusters == clust.to.rm))
  # remove
  DATA <- subset(DATA, !row.names(DATA) %in% clustersToGo)
  attributes(x)$best.clust <- DATA
  # main data
  DATA <- x@main.data
  DATA <- DATA[ , -which(names(DATA) %in% clustersToGo)]
  attributes(x)$main.data <- DATA
  # PCA
  DATA <- x@pca.data
  DATA <- subset(DATA, !row.names(DATA) %in% clustersToGo)
  attributes(x)$pca.data <- DATA
  # tSNE
  DATA <- x@tsne.data
  DATA <- subset(DATA, !row.names(DATA) %in% clustersToGo)
  attributes(x)$tsne.data <- DATA
  #
  DATA <- x@tsne.data.3d
  DATA <- subset(DATA, !row.names(DATA) %in% clustersToGo)
  attributes(x)$tsne.data.3d <- DATA
  return(x)
}
