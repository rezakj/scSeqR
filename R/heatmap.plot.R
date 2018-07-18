#' Create an object of class scSeqR.
#'
#' This function takes data frame and makes an object of class scSeqR.
#' @param x A data frame containing gene counts for cells.
#' @return An object of class scSeqR
#' @examples
#' \dontrun{
#' make.obj(my.data)
#' }
#' @import pheatmap
#' @export
heatmap.plot <- function (x = NULL,
                          gene = "NULL",
                          cluster.by = "clusters",
                          cluster.rows = F,
                          heat.colors = c("blue" ,"white", "red")) {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  ## get main data
  DATAmain <- x@main.data
  AllGenes = row.names(DATAmain)
  absent = which((gene %in% AllGenes) == F)
  absentgenes = gene[absent]
  if(length(absentgenes) != 0)
  {
    absentgenes = paste(absentgenes, collapse=",")
    ToPrint <- paste(absentgenes, "not available in your data.
To see the gene names issue this command: row.names(YOURobject@main.data)", sep=" ")
    stop(print(ToPrint))
  }
  ##### get cluster data
      DATA <- x@best.clust
### get main data
  sub.data <- subset(DATAmain,rownames(DATAmain) %in% gene)
  data.t <- t(sub.data)
  data.expr <- as.data.frame(data.t)
## clusters
  if (cluster.by == "clusters") {
  clusters = DATA
# merge
  mrgd <- merge(clusters, data.expr, by="row.names")
  row.names(mrgd) <- mrgd$Row.names
  mrgd <- mrgd[,-1]
  mrgd <- (mrgd[order(mrgd$clusters, decreasing = F),])
  SideCol <- mrgd[1]
  SideCol$clusters <- sub("^", "cl.",SideCol$clusters)
  data <- mrgd[,-1]
  data <- data.matrix(data)
  data <- log2(data + 1)
  }
  # conditions
  if (cluster.by == "conditions") {
    SideCol <- data.frame(do.call('rbind', strsplit(as.character(rownames(data.expr)),'_',fixed=TRUE)))[1]
    conditions <- as.character(as.matrix(SideCol))
    MyRows = rownames(data.expr)
    conditions <- as.data.frame(cbind(MyRows,conditions = conditions))
    row.names(conditions) <- conditions$MyRows
    conditions <- conditions[2]
    # merge
    mrgd <- merge(conditions, data.expr, by="row.names")
    row.names(mrgd) <- mrgd$Row.names
    mrgd <- mrgd[,-1]
    mrgd <- (mrgd[order(mrgd$conditions, decreasing = F),])
    SideCol <- mrgd[1]
    data <- mrgd[,-1]
    data <- data.matrix(data)
    data <- log2(data + 1)
  }
# Heat map
  mycol <- colorRampPalette(heat.colors)(n = 100)
# return
  return(pheatmap(t(data),
                  col = mycol,
                  show_colnames = F,
                  cluster_rows=cluster.rows,
                  cluster_cols=F,
                  annotation_col = SideCol,
                  scale="row"))
}
