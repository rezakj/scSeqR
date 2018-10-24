setClass("scSeqR", representation (raw.data = "data.frame",
                                   stats = "data.frame",
                                   obj.info = "character",
                                   main.data = "data.frame",
                                   scaled.data = "data.frame",
                                   tsne.data = "data.frame",
                                   umap.data = "data.frame",
                                   pca.data = "data.frame",
                                   pca.info = "list",
                                   opt.pcs = "numeric",
                                   dist.data = "data.frame",
                                   diff.st.data = "data.frame",
                                   tsne.data.3d = "data.frame",
                                   clust.avg = "data.frame",
                                   gene.data = "data.frame",
                                   gene.model = "data.frame",
                                   adt.raw = "data.frame",
                                   adt.main = "data.frame",
                                   clust.cond.freq = "data.frame",
                                   cluster.data = "list",
                                   best.clust = "data.frame",
                                   data.conditions = "character",
                                   norm.factors = "data.frame"))
# hide slots
setMethod("show",
          "scSeqR",
          function(object){
            show(object@obj.info)
          })
