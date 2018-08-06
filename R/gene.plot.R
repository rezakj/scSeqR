#' Find optimal number of clusters.
#'
#' This function takes an object of class scSeqR and finds optimal number of clusters based on three methods.
#' @param x An object of class scSeqR.
#' @param clust.num Number of clusters
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' gene.plot(my.obj, gene = "NULL", box.to.test = 0, box.pval = "sig.signs")
#' }
#' @import ggpubr
#' @export
gene.plot <- function (x = NULL,
                       gene = "NULL",
                       box.to.test = 0,
                       box.pval = "sig.signs",
                       plot.data.type = "tsne",
                       clust.dim = 2,
                       col.by = "clusters",
                       plot.type = "scatterplot",
                       cell.size = 1,
                       cell.colors = c("gray","red"),
                       box.cell.col = "black",
                       box.color = "red",
                       box.line.col = "green",
                       back.col = "white",
                       cell.transparency = 0.5,
                       interactive = TRUE,
                       out.name = "plot") {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  if (gene == "NULL") {
    stop("There is no gene name provided. Please provide a gene name")
  }
  if (length(gene) != 1) {
    stop("currently you can only plot one gene at a time")
  }
  ## get main data
  DATAmain <- x@main.data
  AllGenes = row.names(DATAmain)
  gene.availability = gene %in% AllGenes
  if(gene.availability != TRUE)
  {
    stop("Your gene name is not in the main data. To see the gene names issue this command:
         row.names(YOURobject@main.data)")
  }
  ##### get cluster data
  # 2 dimentions
  if (clust.dim == 2) {
    if (plot.data.type == "tsne") {
      MyTitle = "tSNE Plot"
      DATA <- x@tsne.data
    }
    if (plot.data.type == "pca") {
      MyTitle = "PCA Plot"
      DATA <- x@pca.data
    }
  }
  # 3 dimentions
  if (clust.dim == 3) {
    if (plot.data.type == "tsne") {
      MyTitle = "3D tSNE Plot"
      DATA <- x@tsne.data.3d
    }
    if (plot.data.type == "pca") {
      MyTitle = "3D PCA Plot"
      DATA <- x@pca.data
    }
  }
  # conditions
  if (col.by == "conditions") {
    col.legend <- data.frame(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[1]
    col.legend.box <- factor(as.matrix(col.legend))
  }
  # clusters
  if (col.by == "clusters") {
    if (is.null(x@best.clust)) {
      stop("Clusters are not assigend yet")
    } else {
      col.legend.box <- x@best.clust
      col.legend.box$clusters <- sub("^", "cl.",col.legend.box$clusters)
      col.legend.box <- factor(col.legend.box$clusters)
    }
  }
###### make binary
  # get the gene from the main data
  sub.data <- subset(DATAmain,rownames(DATAmain) == gene)
  data.t <- t(sub.data)
  data.expr <- as.data.frame(data.t)
  data.binary <- data.t > 0
  data.binary <- as.data.frame(data.binary)
  col.legend.bin = data.binary[[gene]]
#### make heamap
  col.legend = log2(data.expr + 1)
  col.legend <- as.numeric(as.matrix(col.legend))
# fix color for ADTs
      if ( length(grep("^ADT_", gene, value = T)) == 1) {
        Lo3=(quantile(col.legend,0.05)) # 1
        Lo2=(quantile(col.legend,0.10)) # 2
        Lo1=(quantile(col.legend,0.25)) # 3
        MID=(quantile(col.legend,0.50)) # 4
        Up1=(quantile(col.legend,0.75)) # 5
        Up2=(quantile(col.legend,0.90)) # 6
        Up3=(quantile(col.legend,0.95)) # 7
        col.legend <- replace(col.legend, col.legend > Up3, Up3)
        col.legend <- replace(col.legend, col.legend > Up2 & col.legend < Up3 ,Up2)
        col.legend <- replace(col.legend, col.legend > Up1 & col.legend < Up2 ,Up1)
        col.legend <- replace(col.legend, col.legend > MID & col.legend < Up1 ,MID)
        col.legend <- replace(col.legend, col.legend > Lo1 & col.legend < MID ,Lo1)
        col.legend <- replace(col.legend, col.legend < Lo2 ,Lo1)
      }
      if ( length(grep("^ADT_", gene, value = T)) == 0) {
        col.legend = col.legend
      }
###
  if (plot.type == "scatterplot") {
  # plot 2d
  if (clust.dim == 2) {
    if (interactive == F) {
      myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                 text = row.names(DATA), color = col.legend)) +
        geom_point(size = cell.size, alpha = cell.transparency) +
        scale_colour_gradient(low = cell.colors[1], high = cell.colors[2], name="") +
        xlab("Dim1") +
        ylab("Dim2") +
        ggtitle(paste(MyTitle,"for (",gene,")")) +
        theme(panel.background = element_rect(fill = back.col, colour = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              legend.key = element_rect(fill = back.col))
    } else {
      myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                 text = row.names(DATA), color = col.legend)) +
        geom_point(size = cell.size, alpha = cell.transparency) +
        scale_colour_gradient(low = cell.colors[1], high = cell.colors[2], name="") +
        xlab("Dim1") +
        ylab("Dim2") +
        ggtitle(paste(MyTitle,"for (",gene,")")) +
        theme_bw()
    }
  }
  # plot 3d
  if (clust.dim == 3) {
    myPLOT <- plot_ly(DATA, x = DATA[,1], y = DATA[,2], z = DATA[,3], text = row.names(DATA),
                      color = col.legend.bin, opacity = cell.transparency, marker = list(size = cell.size + 2)) %>%
      layout(DATA, x = DATA[,1], y = DATA[,2], z = DATA[,3]) %>%
      layout(plot_bgcolor = back.col) %>%
      layout(paper_bgcolor = back.col) %>%
      layout(title = MyTitle,
             scene = list(xaxis = list(title = "Dim1"),
                          yaxis = list(title = "Dim2"),
                          zaxis = list(title = "Dim3")))
    }
  }
  #######
  # plot box plot
  Yaxis1 = data.expr[[gene]]
  Yaxis = Yaxis1 + 1
  Yaxis = log2(Yaxis)
#
  if (plot.type == "boxplot") {
  myPLOT <- ggplot(data.binary, aes(x = col.legend.box, y=Yaxis)) +
    theme_bw() +
    geom_jitter(color = box.cell.col, size = cell.size, alpha = cell.transparency) +
    ggtitle(gene) +
    geom_boxplot(fill = box.color,
                 col = "green",
                 notch = F,
                 outlier.shape = NA,
                 alpha = cell.transparency) +
    ylab("scaled normalized expression") +
    xlab(".")
  # add p-val
  if (box.to.test == 0) {
    AvData <- x@clust.avg
    row.names(AvData) <- AvData$gene
    AvData <- AvData[,-1]
    AvData <- subset(AvData,row.names(AvData) == gene)
    box.to.test <- as.numeric(which.max(AvData))
  }
  if (box.pval == "sig.signs") {
    myPLOT <- myPLOT + stat_compare_means(label = "p.signif", ref.group = box.to.test)
  }
  if (box.pval == "sig.values") {
    myPLOT <- myPLOT + stat_compare_means(aes(label = paste0("p = ", ..p.format..)),ref.group = box.to.test)
  }
  }
  # Bar plot
  if (plot.type == "barplot") {
  myPLOT <- ggplot(data.expr, aes(x = col.legend.box, y = Yaxis1, fill = col.legend.box)) +
    stat_summary(fun.y="mean",
                 geom="bar",
                 alpha = cell.transparency,
                 show.legend = F) +
    ylab("avraged normalized expression") +
    xlab(".") +
    ggtitle(gene) +
    theme_bw()
  }
  # return
  if (interactive == T) {
    OUT.PUT <- paste(out.name, ".html", sep="")
    htmlwidgets::saveWidget(ggplotly(myPLOT), OUT.PUT)
  } else {
    return(myPLOT)
  }
}
