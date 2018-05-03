#' Find optimal number of clusters.
#'
#' This function takes an object of class scSeqR and finds optimal number of clusters based on three methods.
#' @param x An object of class scSeqR.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' pot.clust.num(my.obj)
#' }
#' @export
pot.clust.num <- function (x = NULL) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  require(factoextra)
  require(gridExtra)
  df <- x@tsne.data
  # Elbow method
  Elbow = fviz_nbclust(df, kmeans, method = "wss") + geom_vline(xintercept = 4, linetype = 2) + labs(subtitle = "Elbow method")
  # Silhouette method
  Silhouette = fviz_nbclust(df, kmeans, method = "silhouette") + labs(subtitle = "Silhouette method")
  # Gap statistic
  # nboot = 50 to keep the function speedy.
  # recommended value: nboot= 500 for your analysis.
  # Use verbose = FALSE to hide computing progression.
  set.seed(123)
  Gap.statistic = fviz_nbclust(df, kmeans, nstart = 25,  method = "gap_stat", nboot = 50) + labs(subtitle = "Gap statistic method")
  return(grid.arrange(Elbow, Silhouette,Gap.statistic, nrow = 3))
}
