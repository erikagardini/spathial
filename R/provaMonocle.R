# data<-system.file("extdata", "karlsson-rawcounts.rda", package = "spathial")
# load(data)
#
#
# library(Seurat)
# ### Create Seurat object + Inizialization steps
# seuset<-CreateSeuratObject(counts=rawcounts,min.cells=3,min.features=200)
# normalized_seuset<-NormalizeData(seuset,normalization.method="LogNormalize",scale.factor=10000)
# expmat<-as.matrix(seuset[["RNA"]]@data)
# normalized_expmat <- as.matrix(normalized_seuset[["RNA"]]@data)
# save(expmat,annotation,file="karlsson-expmat.rda")
# save(normalized_expmat,annotation,file="karlsson-normalized_expmat.rda")
#
# expmat<-system.file("extdata", "karlsson-expmat.rda", package = "spathial")
# load(expmat)
#
# cell_annotation <- as.matrix(annotation)
# rownames(cell_annotation) <- colnames(expmat)
# colnames(cell_annotation) <- "cycle"
#
# gene_annotation <- as.matrix(rownames(expmat))
# rownames(gene_annotation) <- rownames(expmat)
# colnames(gene_annotation) <- "gene_short_name"
#
# cds <- new_cell_data_set(expmat,
#                          cell_metadata = cell_annotation,
#                          gene_metadata = gene_annotation)
# cds <- preprocess_cds(cds, num_dim = 100)
# cds <- reduce_dimension(cds, preprocess_method = "PCA")
# cds <- cluster_cells(cds)
# cds <- learn_graph(cds)
# png("monocle3_graph.png",w=3000,h=3000,res=400)
# plot_cells(cds,
#            color_cells_by = "cycle",
#            label_groups_by_cluster=FALSE,
#            label_leaves=FALSE,
#            label_branch_points=FALSE)
# dev.off()
# ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
# save(ciliated_cds_pr_test_res, file="monocle_res.rda")
#
# normalized_expmat<-system.file("extdata", "karlsson-normalized_expmat.rda", package = "spathial")
# load(normalized_expmat)
# X <- t(normalized_expmat)
# X_labels <- annotation
#
# boundary_init <- spathial::spathialBoundaryIds(X, X_labels, mode=2, from="G1", to="G2/M")
# spathial_res <- spathial::spathialWay(boundary_init$X, boundary_init$boundary_ids, NC = 50)
# spathial::spathialPlot(boundary_init$X, boundary_init$X_labels, boundary_init$boundary_ids, spathial_res, perplexity = 30, mask = NULL)
# stat_spathial <- spathial::spathialStatistics(spathial_res)
#
# spathial_best_genes <- names(stat_spathial$correlations[which(stat_spathial$correlations > 0.8 | stat_spathial$correlations < -0.8)])
# monocle_best_genes <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
#
# genes<-system.file("extdata", "karlsson-genes.rda", package = "spathial")
# load(genes)
# common<-intersect(spathial_best_genes, monocle_best_genes)
# best_spathial <- setdiff(spathial_best_genes, monocle_best_genes)
# best_monocle <- setdiff(monocle_best_genes, spathial_best_genes)
# best_spathial[which(best_spathial %in% group_2)]
# best_spathial[which(best_spathial %in% group_3)]
#
