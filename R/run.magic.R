#' Run UMAP on PCA data
#'
#' This function takes an object of class scSeqR and runs UMAP on PCA data.
#' @param x An object of class scSeqR.
#' @param dims PC dimentions to be used for UMAP analysis.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' my.obj <- run.umap(my.obj, dims = 1:10)
#' }
#' @import umap
#' @export
run.magic <- function (x = NULL) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
### load packages
  # https://github.com/KrishnaswamyLab/MAGIC/tree/master/Rmagic
  require(Rmagic)
  require(viridis)
  require(phateR)
  # get data
  DATA <- x@main.data
  DATA <- as.data.frame(t(DATA))
  data_MAGIC <- magic(DATA, genes="all_genes")
  DATA <- as.data.frame(t(data_MAGIC$result))
  # return
  attributes(x)$main.data <- DATA
  # return
  return(x)
}
