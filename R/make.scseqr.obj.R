#' Create an object of class scSeqR.
#'
#' This function takes data frame and makes an object of class scSeqR.
#' @param x A data frame containing gene counts for cells.
#' @return An object of class scSeqR
#' @examples
#' \dontrun{
#' make.obj(my.data)
#' }
#' @export
make.scseqr.obj <- function (x = NULL) {
  setClass("scSeqR", representation (raw.data = "data.frame",
                                     row.names = "character",
                                     stats = "data.frame",
                                     main.data = "data.frame",
                                     scaled.data = "data.frame",
                                     tsne.data = "data.frame"))
  myOBJ <- new("scSeqR", raw.data = x)
  attributes(myOBJ)$row.names <- row.names(x)
  return(myOBJ)
}
