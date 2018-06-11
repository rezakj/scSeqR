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
find.marker <- function (x = NULL,
                         fold.change = 1.5) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  # load data
  DATA <- x@clust.avg
  # for each column divide by rest and make FC and anything less then FC put NA
  A = c(1,4,7,5)

  for (i in A) {
  FC <- i / mean(A[-which((A) %in% i)])
  log2FC <- log2(FC)
  }

  return(object)
}
