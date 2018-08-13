#' Find optimal number of PCs for clustering
#'
#' This function takes an object of class scSeqR and finds optimal number of PCs for clustering.
#' @param x An object of class scSeqR.
#' @param pcs.in.plot Number of PCs to show in plot, defult = 50.
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' find.opt.pcs(my.obj)
#' }
#'
#' @export
find.opt.pcs <- function (x = NULL, pcs.in.plot = 50) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  DATA <- x@pca.info$sdev
  OPTpcs <- mean(DATA)*2
  OPTpcs <- (DATA > OPTpcs)
  OPTpcs <- length(OPTpcs[OPTpcs==TRUE])
  Titel <- paste("Optimal number of PCs (1:", OPTpcs, ")", sep="")
  # fix DATA
  DATA <- as.data.frame(DATA)
  w = log1p(1:100)
  Rows <- c(1:dim(DATA)[1])
  DATA <- cbind(Rows, DATA)
  colnames(DATA) <- c("PCs","SDs")
  DATA <- head(DATA,pcs.in.plot)
### plot
  myPLOT <- ggplot(DATA,aes(y=SDs,x=PCs)) +
    geom_point() + geom_line() +
    scale_x_continuous(trans = "log1p") +
    scale_y_continuous(trans = "log1p") +
    geom_vline(xintercept = OPTpcs,linetype="dotted")+
    theme_bw(base_size = 16)+
    xlab("Principal Components") +
    ylab("PC Standard Deviations") +
    ggtitle(Titel)
# return
  return(myPLOT)
}

