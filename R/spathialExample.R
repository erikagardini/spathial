# Set Random Seed
set.seed(1)
# Number of optimization variables i.e. path waypoints
NC <- 50

myfile<-system.file("extdata", "liver_tcga.csv", package = "spathial")
data<-read.csv(myfile,as.is=TRUE,header=TRUE)
X <- data[,2:(dim(data)[2]-1)]
X_labels <- data[,dim(data)[2]]
rownames(X)<-data[,1]
rownames(X)<-gsub("\\.","-",rownames(X))

# # SUBSET OF ENTRIES, JUST FOR TESTING
# X <- X[430:450,1:5]
# X_labels <- X_labels[430:450]

#Variance filtering
varFiltering <- TRUE
nvars <- 1000
if(varFiltering){
  vars <- apply(X,2,var)
  topfeatures <- names(sort(vars, dec=TRUE))[1:nvars]
  X <- X[,topfeatures]
}

# Choose the starting and the ending points
boundary_init <- spathialBoundaryIds(X, X_labels, mode=1, from=2, to=1)
boundary_ids <- boundary_init$boundary_ids
X <- boundary_init$X
X_labels <- boundary_init$X_labels

# Prefilter data
filter_res <- spathialPrefiltering(X, X_labels, boundary_ids)

X_filtered <- filter_res$X_filtered
X_labels_filtered <- filter_res$X_labels_filtered
X_garbage <- filter_res$X_garbage
X_labels_garbage <- filter_res$X_labels_garbage
boundary_ids_filtered <- filter_res$boundary_ids

# X_filtered <- X
# X_labels_filtered <- X_labels
# boundary_ids_filtered <- boundary_ids

# Compute spathial
spathial_res <- spathialWay(X_filtered, X_labels_filtered, boundary_ids_filtered, NC, neighbors = 1)

#Labels for each waypoint with knn
ppath_labels <- spathialLabels(X_filtered, X_labels_filtered, spathial_res)
#Plot the path in 2D using Rtsne
#spathialPlot2D(X_filtered, X_labels_filtered, boundary_ids, spathial_res, perplexity_value=30, X_garbage, X_labels_garbage)
spathialPlot2D(X_filtered, X_labels_filtered, boundary_ids, spathial_res, perplexity_value=30, X_garbage, X_labels_garbage)

# save(spathial_res, file="res_parallelized.rda")
#
# load(file = "~/spathial/spathial_res.rda")
#
# perturbed_path <- spathial_res
# principal_path <- spathial_res
# principal_path$perturbed_path <- NULL
#
# #Correlation along the path
# message("SPATHIAL CORRELATION")
# corr_values_perturbed_path <- spathial_corr(perturbed_path)
# corr_perturbed_path <- unlist(corr_values_perturbed_path, use.names=TRUE)
#
# corr_values_principal_path <- spathial_corr(principal_path)
# corr_principal_path <- unlist(corr_values_principal_path, use.names=TRUE)
#
# library(msigdbr)
# library(fgsea)
# library(ggplot2)
# m_df = msigdbr(species = "Homo sapiens", category = "H")
# genes <- m_df %>% split(x=.$gene_symbol, f=.$gs_name)
#
# fgseaRes_perturbed <- fgsea(pathways = genes,
#                   stats = corr_perturbed_path,
#                   minSize=5,
#                   maxSize=500,
#                   nperm=100000)
#
# plotEnrichment(genes[["HALLMARK_P53_PATHWAY"]],
#                corr_perturbed_path) + labs(title="HALLMARK_P53_PATHWAY")
#
# fgseaRes_principal <- fgsea(pathways = genes,
#                   stats = corr_principal_path,
#                   minSize=5,
#                   maxSize=500,
#                   nperm=100000)
#
# plotEnrichment(genes[["HALLMARK_P53_PATHWAY"]],
#                corr_principal_path) + labs(title="HALLMARK_P53_PATHWAY")
