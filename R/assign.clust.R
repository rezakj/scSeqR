#' Find optimal number of clusters.
#'
#' This function takes an object of class scSeqR and finds optimal number of clusters based on three methods.
#' @param x An object of class scSeqR.
#' @param clust.num Number of clusters
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' assign.clust(my.obj, clust.num = 7)
#' }
#' @export
assign.clust <- function (x = NULL, clust.num = 0) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  if ( clust.num == 0) {
    stop("please provite the optimal number of clusters")
  }
  df <- x@tsne.data
  DATA <- x@tsne.data
  row.names(df) <- df$clls
  df <- df[,-1]
  d_tsne_1=df
  d_tsne_1_original <- d_tsne_1
  ClustNum = clust.num
  ## Creating k-means clustering model, and assigning the result to the data used to create the tsne
  fit_cluster_kmeans=kmeans(scale(d_tsne_1), ClustNum)
  d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)
  ## Creating hierarchical cluster model, and assigning the result to the data used to create the tsne
  fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1)))
  ## setting 3 clusters as output
  d_tsne_1_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=ClustNum))
  ###
  tsne.data <- cbind(cells = DATA$cells,d_tsne_1_original)
  attributes(x)$tsne.data <- tsne.data
  return(x)
  }
