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
make.obj <- function (x = NULL) {
  setClass("scSeqR", representation (raw.data = "data.frame",
                                     row.names = "character",
                                     stats = "data.frame",
                                     main.data = "data.frame",
                                     scaled.data = "data.frame",
                                     tsne.data = "data.frame"))
  object <- new(Class = "scSeqR",raw.data = x)
  return(object)
}
