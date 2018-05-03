#' Scale data
#'
#' This function takes an object of class scSeqR and scales the normalized data.
#' @param x An object of class scSeqR.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' scale.data(my.obj)
#' }
#' @export
scale.data <- function (x = NULL) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  DATA <- x@main.data
  NormLog <- log(DATA+0.1)
  NormLogScale <- as.data.frame(scale(NormLog))
  attributes(x)$scaled.data <- NormLogScale
  return(x)
}

