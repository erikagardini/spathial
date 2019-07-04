# Set Random Seed
set.seed(1)
# Number of optimization variables i.e. path waypoints
NC <- 50

myfile<-system.file("extdata", "thyroid_tcga.csv", package = "spathial")
data<-read.csv(myfile,as.is=TRUE,header=TRUE)
X <- data[,2:(dim(data)[2]-1)]
X_labels <- data[,dim(data)[2]]
rownames(X)<-data[,1]
rownames(X)<-gsub("\\.","-",rownames(X))

#SUBSET OF ENTRIES, JUST FOR TESTING
#X <- X[438:445,1:5]
#X_labels <- X_labels[438:445]

#Variance filtering
# varFiltering <- TRUE
# nvars <- 10000
# if(varFiltering){
#   vars <- apply(X,2,var)
#   topfeatures <- names(sort(vars, dec=TRUE))[1:nvars]
#   X <- X[,topfeatures]
# }

#Choose the starting and the ending points
boundary_init <- spathial_boundary_ids(X, X_labels, mode=1, from=2, to=1)
boundary_ids <- boundary_init$boundary_ids
X <- boundary_init$X
X_labels <- boundary_init$X_labels

#Compute spathial
spathial_res <- spathial_way_multiple(X, X_labels, boundary_ids, NC, prefiltering=FALSE, negb = 2)
#Labels for each waypoint with knn
ppath_labels <- spathial_labels(X, X_labels, spathial_res)
#Plot the path in 2D using tsne
spathial_2D_plot(X, X_labels, boundary_ids, spathial_res, 30)

save(spathial_res, file="res_parallelized.rda")

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
