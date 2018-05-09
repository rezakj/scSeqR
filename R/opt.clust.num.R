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
pot.clust.num <- function (x = NULL, max.clust = 15, gap.stat.nboot = 100, verbose = TRUE) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  require(factoextra)
  require(gridExtra)
  df <- x@tsne.data
  row.names(df) <- df$cells
  df <- df[,-1]
  df <- df[1:2]
  # Elbow method
  Elbow = fviz_nbclust(df, kmeans, method = "wss", k.max = max.clust) +
    geom_vline(xintercept = 4, linetype = 2) +
    labs(subtitle = "Elbow method")
  # Silhouette method
  Silhouette = fviz_nbclust(df, kmeans, method = "silhouette", k.max = max.clust) +
    labs(subtitle = "Silhouette method")
  # Gap statistic
  # nboot = 50 to keep the function speedy.
  # recommended value: nboot= 500 for your analysis.
  # Use verbose = FALSE to hide computing progression.
  set.seed(123)
  Gap.statistic = fviz_nbclust(df, kmeans,
                               nstart = 25,
                               method = "gap_stat",
                               nboot = gap.stat.nboot,
                               k.max = max.clust) +
    labs(subtitle = "Gap statistic method")
  return(grid.arrange(Elbow, Silhouette,Gap.statistic, nrow = 3))
}
