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
  INFO = "An object of class scSeqR version:"
  INFO = paste(INFO, packageVersion("scSeqR"))
  Data.Dim = dim(x)
  Data.Dim <- paste(Data.Dim , collapse=",")
  Data.Dim <- paste("Raw/original data dimentions (rows,columns):", Data.Dim)
  DATA <- colnames(x)
  Col.n <- head(DATA,3)
  Col.n <- paste(Col.n, collapse=",")
  Col.n <- paste("Columns names:" , Col.n, "...")
  Row.n <- head(row.names(x),3)
  Row.n <- paste(Row.n, collapse=",")
  Row.n <- paste("Row names:" , Row.n, "...")
  # get conditions
  do <- data.frame(do.call('rbind', strsplit(as.character(head(DATA,1)),'_',fixed=TRUE)))
  do <- dim(do)[2]
  if (do == 2) {
    My.Conds <- data.frame(do.call('rbind', strsplit(as.character(DATA),'_',fixed=TRUE)))[1]
    My.Conds <- as.character(unique(as.matrix(My.Conds)))
    My.Conds <- paste(My.Conds , collapse=",")
    My.Conds <- paste("Data conditions:", My.Conds)
  } else {
    My.Conds = "Data conditions: no conditions/single sample"
  }
# paste
INFO.to.show <- paste(INFO, Data.Dim, My.Conds, Col.n, Row.n, sep="\n")
INFO.to.show <- capture.output(cat(INFO.to.show))
# make object
  object <- new(Class = "scSeqR", obj.info = INFO.to.show, raw.data = x, data.conditions = My.Conds)
# return
  return(object)
}
