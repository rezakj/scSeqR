#' Run tSNE on PCA data
#'
#' This function takes an object of class scSeqR and runs tSNE on PCA data.
#' @param x An object of class scSeqR.
#' @param dims PC dimentions to be used for tSNE analysis.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' my.obj <- run.pc.tsne(my.obj, dims = 1:10)
#' }
#' @import umap
#' @export
run.umap <- function (x = NULL,
                         dims = 1:20) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  # get PCA data
  DATA <- x@pca.data
  TransPosed <- DATA[dims]
  #  2 dimention
  #  if (clust.dim == 2) {
  # TransPosed <- t(TopNormLogScale)
  myUMAP = umap(TransPosed)
  myUMAP = as.data.frame((myUMAP$layout))
  attributes(x)$umap.data <- myUMAP
# return
  return(x)
}
