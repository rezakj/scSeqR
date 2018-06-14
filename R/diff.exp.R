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
diff.exp <- function (x = NULL,
                      plot.type = "tsne",
                      clust.dim = 2,
                      de.by = "clusters",
                      cond.1 = "array",
                      cond.2 = "array") {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  if (clust.dim != 2 && clust.dim != 3) {
    stop("clust.dim should be either 2 or 3")
  }
  ###########
  dat <- x@main.data
  # 2 dimentions
  if (clust.dim == 2) {
    if (plot.type == "tsne") {
      MyTitle = "tSNE Plot"
      DATA <- x@tsne.data
    }
    if (plot.type == "pca") {
      MyTitle = "PCA Plot"
      DATA <- x@pca.data
    }
  }
  # 3 dimentions
  if (clust.dim == 3) {
    if (plot.type == "tsne") {
      MyTitle = "3D tSNE Plot"
      DATA <- x@tsne.data.3d
    }
    if (plot.type == "pca") {
      MyTitle = "3D PCA Plot"
      DATA <- x@pca.data.3d
    }
  }
  ############## set wich clusters you want as condition 1 and 2
  CondA = cond.1
  CondB = cond.2
  CondAnames = paste(CondA, collapse="_")
  CondBnames = paste(CondB, collapse="_")
  # conditions
  if (de.by == "conditions") {
    col.legend <- data.frame(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[1]
    conditions <- as.character(as.matrix(col.legend))
    IDs = as.character(rownames(DATA))
    myconds <- as.data.frame(cbind(IDs = IDs, conditions = conditions))
    rownames(myconds) <- myconds$IDs
    Table <- myconds[2]
    Cluster0 <- row.names(subset(Table, Table$conditions %in% CondA))
    Cluster1 <- row.names(subset(Table, Table$conditions %in% CondB))
  }
  # clusters
  if (de.by == "clusters") {
    if (is.null(DATA$clusters)) {
      stop("Clusters are not assigend yet, please run assign.clust fisrt.")
    } else {
      Table=DATA
      Cluster0 <- row.names(subset(Table, Table$clusters %in% CondA))
      Cluster1 <- row.names(subset(Table, Table$clusters %in% CondB))
    }
  }
  #
  ############## Filter
  cond1 <- dat[,Cluster0]
  cond2 <- dat[,Cluster1]
  ### merge both for pval length not matching error
  mrgd <- merge(cond1, cond2, by="row.names")
  row.names(mrgd) <- mrgd$Row.names
  mrgd <- mrgd[,-1]
  mrgd <- data.matrix(mrgd)
  # mean
  meansCond1 <- apply(cond1, 1, function(cond1) {mean(cond1)})
  meansCond2 <- apply(cond2, 1, function(cond2) {mean(cond2)})
  baseMean <- apply(mrgd, 1, function(mrgd) {mean(mrgd)})
  # FC
  FC <- meansCond2/meansCond1
  FC.log2 <- log2(FC)
  # dims
  Cond1_Start <- 1
  Cond1_End <- dim(cond1)[2]
  Cond2_Start <- dim(cond1)[2] + 1
  Cond2_End <- dim(cond1)[2] + dim(cond2)[2]
  # pval
  Pval <- apply(mrgd, 1, function(mrgd) {
    t.test(x = mrgd[Cond1_Start:Cond1_End], y = mrgd[Cond2_Start:Cond2_End])$p.value
  })
  # padj
  FDR <- p.adjust(Pval)
  # combine
  Stats <- cbind(
    baseMean = baseMean,
    MeanCond1 = meansCond1,
    MeanCond2 = meansCond2,
    FC=FC,
    FC.log2=FC.log2,
    t.test=Pval,
    t.test.adj=FDR)
  # name column
  colnames(Stats) <- c("baseMean",CondAnames, CondBnames,"foldChange","log2FoldChange","pval","padj")
# return
    return(Stats)
}

