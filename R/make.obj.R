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
  # get info
  INFO = "Class: an object of class scSeqR"
  Data.Dim = dim(x)
  DATA <- colnames(x)
  # get conditions
  do <- data.frame(do.call('rbind', strsplit(as.character(head(DATA,1)),'_',fixed=TRUE)))
  do <- dim(do)[2]
  if (do == 2) {
    My.Conds <- data.frame(do.call('rbind', strsplit(as.character(DATA),'_',fixed=TRUE)))[1]
    My.Conds <- as.character(unique(as.matrix(My.Conds)))
  } else {
    My.Conds = "No conditions or single sample"
  }
# paste
Data.Dim <- paste(Data.Dim , collapse=",")
Data.Dim <- paste("Raw/original data dimentions (rows,columns):", Data.Dim)
My.Conds <- paste(My.Conds , collapse=",")
My.Conds <- paste("Data conditions:", My.Conds)
INFO.to.show <- paste(INFO, Data.Dim, My.Conds, sep="\n")
INFO.to.show <- capture.output(cat(INFO.to.show))
# make object
  object <- new(Class = "scSeqR", obj.info = INFO.to.show, raw.data = x)
# return
  return(object)
}
