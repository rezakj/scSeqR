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
add.adt <- function (x = NULL, adt.data = "data.frame") {
  if (class(adt.data) != "data.frame") {
    stop("ADT data should be a data frame object")
  }
  attributes(x)$adt.raw <- adt.data
  return(x)
}
