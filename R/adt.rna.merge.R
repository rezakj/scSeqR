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
adt.rna.merge <- function (x = NULL, adt.data = "raw") {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  ### get ADT data
  if (adt.data == "raw") {
  DATA <- as.data.frame(t(x@adt.raw))
  }
  if (adt.data == "main") {
  DATA <- as.data.frame(t(x@adt.main))
  }
  ### get RNA data
  DATArna <- as.data.frame(t(x@main.data))
  ### merge
  merged.data <- merge(DATA, DATArna, by="row.names", all.x=F, all.y=T)
  rownames(merged.data) <- merged.data$Row.names
  merged.data <- merged.data[,-1]
  merged.data <- as.data.frame(t(merged.data))
  attributes(x)$main.data <- merged.data
  return(x)
}
