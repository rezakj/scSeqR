#' Cluster based on tSNE
#'
#' This function takes an object of class scSeqR and clusters the normalized data.
#' @param x An object of class scSeqR.
#' @param clust.method A method to choose your genes for clustering. Choose from "base.mean.rank" or "gene.model", defult is "base.mean.rank".
#' @param top.rank Number of ranked genes you want to use for clustering, defult = 500.
#' @param clust.type Choose from "tsne", "pca" or "distance".
#' @param dist.method Choose "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski", defult = "euclidean".
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' cluster.data(my.obj, clust.method = "base.mean.rank", top.rank = 500, clust.dim = 3)
#' }
#' @import Rtsne
#' @export
run.tsne <- function (x = NULL,
                          clust.method = "base.mean.rank",
                          top.rank = 500,
                          gene.list = "my_model_genes.txt") {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
#  if (clust.dim != 2 && clust.dim != 3) {
#    stop("clust.dim should be either 2 or 3")
#  }
  if (clust.method == "dispersed.genes" && clust.method == "both") {
    stop("dispersed.genes and both are not implemented yet")
  }
  # geth the genes and scale them based on model
  DATA <- x@main.data
  # model base mean rank
  if (clust.method == "base.mean.rank") {
    raw.data.order <- DATA[ order(rowMeans(DATA), decreasing = T), ]
    topGenes <- head(raw.data.order,top.rank)
    NormLog <- log(topGenes+0.1)
    TopNormLogScale <- as.data.frame(scale(NormLog))
  }
  # gene model
  if (clust.method == "gene.model") {
    if (!file.exists(gene.list)) {
      stop("please provide a file with gene names for clustering")
    } else {
      genesForClustering <- readLines(gene.list)
      topGenes <- subset(DATA, rownames(DATA) %in% genesForClustering)
      NormLog <- log(topGenes+0.1)
      TopNormLogScale <- as.data.frame(scale(NormLog))
    }
  }
#  2 dimention
#  if (clust.dim == 2) {
      TransPosed <- t(TopNormLogScale)
      tsne <- Rtsne(TransPosed, check_duplicates = FALSE, dims = 2)
      tsne.data = as.data.frame(tsne$Y)
      tsne.data = cbind(cells = row.names(TransPosed),tsne.data)
      rownames(tsne.data) <- tsne.data$cells
      tsne.data <- tsne.data[,-1]
      attributes(x)$tsne.data <- tsne.data
#  }
# choose 3 demention
  # tSNE
      TransPosed <- t(TopNormLogScale)
      tsne <- Rtsne(TransPosed, check_duplicates = FALSE, dims = 3)
      tsne.data = as.data.frame(tsne$Y)
      tsne.data = cbind(cells = row.names(TransPosed),tsne.data)
      rownames(tsne.data) <- tsne.data$cells
      tsne.data <- tsne.data[,-1]
      attributes(x)$tsne.data.3d <- tsne.data
  return(x)
}
