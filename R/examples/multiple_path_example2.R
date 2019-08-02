# ####Load input data
# myfile<-system.file("extdata", "thyroid_tcga.csv", package = "spathial")
# data<-read.csv(myfile,as.is=TRUE,header=TRUE)
#
# # #Variance filtering
#
# X <- data[,2:(ncol(data)-1)]
# rownames(X)<-paste0("sam",rownames(X))
# X_labels <- data$Y
#
# #### Choose the starting and the ending points
# boundary_init <- spathialBoundaryIds(X, X_labels, mode=2, from=2, to=1)
# boundary_ids <- boundary_init$boundary_ids
# X <- boundary_init$X
# X_labels <- boundary_init$X_labels
#
# # Compute spathial
# NC <- 50
# spathial_res_perturbed <- spathialWay(X, X_labels, boundary_ids, NC, neighbors = 2)
# spathial_res <- spathialWay(X, X_labels, boundary_ids, NC, neighbors = 0)
#
# #Plot the path in 2D using Rtsne
# spathialPlot(X, X_labels, boundary_ids, spathial_res, perplexity_value=30)
# spathialPlot(X, X_labels, boundary_ids, spathial_res_perturbed, perplexity_value=30)
#
# # ####Correlation along the path
# statistics_spathial_res <- spathialStatistics(spathial_res)
# statistics_spathial_res_perturbed <- spathialStatistics(spathial_res_perturbed)
#
# pp_corr <- statistics_spathial_res$correlations
# pp_corr <- unlist(pp_corr, use.names=TRUE)
#
# ppp_corr <- statistics_spathial_res_perturbed$correlations
# ppp_corr <- unlist(ppp_corr, use.names=TRUE)
#
# library(msigdbr)
# library(fgsea)
# library(ggplot2)
# m_df = msigdbr(species = "Homo sapiens", category = "H")
# genes <- m_df %>% split(x=.$gene_symbol, f=.$gs_name)
#
# # Custom subset of pathways
# m_all = msigdbr(species = "Homo sapiens")
# genes_all <- m_all %>% split(x=.$gene_symbol, f=.$gs_name)
# #genes_thyroid <- genes_all[grep("THYROID",names(genes_all))]
#
# fgseaRes_pp <- fgsea(pathways = genes_all,
#                      stats = pp_corr,
#                      minSize=5,
#                      maxSize=500,
#                      nperm=100000)
#
# plotEnrichment(genes[["HALLMARK_P53_PATHWAY"]],
#                pp_corr) + labs(title="HALLMARK_P53_PATHWAY")
#
# fgseaRes_ppp <- fgsea(pathways = genes_all,
#                       stats = ppp_corr,
#                       minSize=5,
#                       maxSize=500,
#                       nperm=100000)
#
# plotEnrichment(genes[["HALLMARK_P53_PATHWAY"]],
#                ppp_corr) + labs(title="HALLMARK_P53_PATHWAY")
