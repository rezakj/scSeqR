#' Plot nGenes, UMIs and perecent mito
#'
#' This function takes an object of class scSeqR and creats plot.
#' @param x An object of class scSeqR.
#' @param cell.size A numeric value for the size of the cells, defults = 1.
#' @param plot.type Choose between "tsne" and "pca", defult = "tsne".
#' @param cell.color Choose cell color if col.by = "monochrome", defult = "black".
#' @param back.col Choose background color, defult = "black".
#' @param col.by Choose between "clusters", "conditions" or "monochrome", defult = "clusters".
#' @param cell.transparency A numeric value between 0 to 1, defult = 0.5.
#' @param clust.dim A numeric value for plot dimentions. Choose either 2 or 3, defult = 2.
#' @param interactive If TRUE an html intractive file will be made, defult = TRUE.
#' @param out.name Output name for html file if interactive = TRUE defult = "plot".
#' @return An object of class scSeqR.
#' @examples
#' \dontrun{
#' tsne.plot(my.obj)
#' }
#' @import ggplot2
#' @export
cluster.plot <- function (x = NULL,
                          cell.size = 1,
                          plot.type = "tsne",
                          cell.color = "black",
                          back.col = "white",
                          col.by = "clusters",
                          cell.transparency = 0.5,
                          clust.dim = 2,
                          density = F,
                          interactive = TRUE,
                          out.name = "plot") {
  if ("scSeqR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  if (clust.dim != 2 && clust.dim != 3) {
    stop("clust.dim should be either 2 or 3")
  }
  ###########
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
      DATA <- x@pca.data
    }
  }
  # conditions
  if (col.by == "conditions") {
    col.legend <- data.frame(do.call('rbind', strsplit(as.character(rownames(DATA)),'_',fixed=TRUE)))[1]
    col.legend <- factor(as.matrix(col.legend))
    }
  # clusters
  # always use hierarchical (k means changes everytime you run)
  if (col.by == "clusters") {
    if (is.null(x@cluster.data$Best.partition)) {
      stop("Clusters are not assigend yet, please run assign.clust fisrt.")
    } else {
      col.legend <- factor(x@best.clust$clusters)
    }
  }
  # monochrome
  if (col.by == "monochrome") {
    if (clust.dim == 2) {
    myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                               text = row.names(DATA))) +
      geom_point(size = cell.size, col = cell.color, alpha = cell.transparency) +
      xlab("Dim1") +
      ylab("Dim2") +
      ggtitle(MyTitle) +
      scale_color_discrete(name="") +
      theme(panel.background = element_rect(fill = back.col, colour = "black"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            legend.key = element_rect(fill = back.col)) + theme_bw()
    }
    if (clust.dim == 3) {
      myPLOT <- plot_ly(DATA, x = DATA[,1], y = DATA[,2], z = DATA[,3], text = row.names(DATA),
              opacity = cell.transparency, marker = list(size = cell.size + 2, color = cell.color)) %>%
        layout(DATA, x = DATA[,1], y = DATA[,2], z = DATA[,3]) %>%
        layout(plot_bgcolor = back.col) %>%
        layout(paper_bgcolor = back.col) %>%
        layout(title = MyTitle,
             scene = list(xaxis = list(title = "Dim1"),
                          yaxis = list(title = "Dim2"),
                          zaxis = list(title = "Dim3")))
      }
    }
# plot 2d
  if (clust.dim == 2) {
    if (interactive == F) {
      myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                 text = row.names(DATA), color = col.legend)) +
        geom_point(size = cell.size, alpha = cell.transparency) +
        xlab("Dim1") +
        ylab("Dim2") +
        ggtitle(MyTitle) +
        scale_color_discrete(name="") +
        theme(panel.background = element_rect(fill = back.col, colour = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              legend.key = element_rect(fill = back.col))
    } else {
      myPLOT <- ggplot(DATA, aes(DATA[,1], y = DATA[,2],
                                 text = row.names(DATA), color = col.legend)) +
        geom_point(size = cell.size, alpha = cell.transparency) +
        xlab("Dim1") +
        ylab("Dim2") +
        ggtitle(MyTitle) +
        scale_color_discrete(name="") +
        theme_bw()
    }
  }
# plot 3d
  if (clust.dim == 3) {
    myPLOT <- plot_ly(DATA, x = DATA[,1], y = DATA[,2], z = DATA[,3], text = row.names(DATA),
                      color = col.legend, opacity = cell.transparency, marker = list(size = cell.size + 2)) %>%
      layout(DATA, x = DATA[,1], y = DATA[,2], z = DATA[,3]) %>%
      layout(plot_bgcolor = back.col) %>%
      layout(paper_bgcolor = back.col) %>%
      layout(title = MyTitle,
           scene = list(xaxis = list(title = "Dim1"),
                        yaxis = list(title = "Dim2"),
                        zaxis = list(title = "Dim3")))
  }
  # density plot
  if (density == T) {
    myPLOT <- ggplot(DATA, aes(DATA[,2], fill = col.legend)) +
      geom_density(alpha=cell.transparency) +
      xlab("Dim2") +
      scale_color_discrete(name="") + theme_bw()
  }
# return
  if (interactive == T) {
    OUT.PUT <- paste(out.name, ".html", sep="")
    htmlwidgets::saveWidget(ggplotly(myPLOT), OUT.PUT)
  } else {
    return(myPLOT)
  }
}

