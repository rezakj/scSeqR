#' Scale data
#'
#' This function takes an object of class scSeqR and scales the normalized data.
#' @param x An object of class scSeqR.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' data.scale(my.obj)
#' }
#'
#' @export
data.scale <- function (x = NULL) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  DATA <- x@main.data
  DATA1 = DATA + 1
  NormLog = log2(DATA1)
  attributes(x)$scaled.data <- NormLog
  return(x)
}

