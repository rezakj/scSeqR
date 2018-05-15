#' Cluster based on tSNE
#'
#' This function takes an object of class scSeqR and clusters the normalized data.
#' @param x An object of class scSeqR.
#' @param clust.method A method to choose your genes for clustering. There are three options "base.mean.rank","dispersed.genes","both" or "my.genes",defults is "base.mean.rank".
#' @param top.rank Number of ranked genes you want to use for clustering, defult = 500.
#' @param clust.dim Number of dimentions you want to consider for clustering, defult =2.
#' @param clust.type Choose from "tsne", "pca" or "distance".
#' @param dist.method Choose "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski", defult = "euclidean".
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' cluster.data(my.obj, clust.method = "base.mean.rank", top.rank = 500, clust.dim = 3)
#' }
#' @import Rtsne
#' @export
cluster.data <- function (x = NULL,
                          clust.method = "base.mean.rank",
                          top.rank = 500,
                          gene.list = NULL,
                          clust.dim = 2,
                          clust.type = "tsne",
                          dist.method = "euclidean") {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  if (clust.dim != 2 && clust.dim != 3) {
    stop("clust.dim should be either 2 or 3")
  }
  if (clust.method == "dispersed.genes" && clust.method == "both") {
    stop("dispersed.genes and both are not implemented yet")
  }
  DATA <- x@main.data
  if (clust.method == "base.mean.rank") {
    raw.data.order <- DATA[ order(rowMeans(DATA), decreasing = T), ]
    topGenes = head(raw.data.order,top.rank)
    NormLog <- log(topGenes+0.1)
    TopNormLogScale <- as.data.frame(scale(NormLog))
  }
  if (clust.dim == 2) {
    if (clust.type == "tsne") {
      TransPosed <- t(TopNormLogScale)
      tsne <- Rtsne(TransPosed, check_duplicates = FALSE, dims = clust.dim)
      tsne.data = as.data.frame(tsne$Y)
      tsne.data = cbind(cells = row.names(TransPosed),tsne.data)
      rownames(tsne.data) <- tsne.data$cells
      tsne.data <- tsne.data[,-1]
      attributes(x)$tsne.data <- tsne.data
    }
    if (clust.type == "pca") {
      counts.pca <- prcomp(TopNormLogScale,center = T,scale. = T)
      dataPCA = data.frame(counts.pca$rotation)[1:2]
      attributes(x)$pca.data <- dataPCA
    }
  }
  if (clust.dim == 3) {
    if (clust.type == "tsne") {
      TransPosed <- t(TopNormLogScale)
      tsne <- Rtsne(TransPosed, check_duplicates = FALSE, dims = clust.dim)
      tsne.data = as.data.frame(tsne$Y)
      tsne.data = cbind(cells = row.names(TransPosed),tsne.data)
      rownames(tsne.data) <- tsne.data$cells
      tsne.data <- tsne.data[,-1]
      attributes(x)$tsne.data.3d <- tsne.data
    }
    if (clust.type == "pca") {
      counts.pca <- prcomp(TopNormLogScale,center = T,scale. = T)
      dataPCA = data.frame(counts.pca$rotation)[1:3]
      attributes(x)$pca.data.3d <- dataPCA
    }
  }
    if (clust.type == "distance") {
  dists = dist(t(TopNormLogScale), method = dist.method)
  dists.data = as.data.frame(as.matrix(dists))
  attributes(x)$dist.data <- dists.data
  }
  return(x)
}

