# Set Random Seed
set.seed(1)
# Number of optimization variables i.e. path waypoints
NC <- 50

####Init input data
myfile<-system.file("extdata", "liver_tcga.csv", package = "spathial")
data<-read.csv(myfile,as.is=TRUE,header=TRUE)
data$Y[which(data$Y == 'G1')] <- 1
data$Y[which(data$Y == 'S')] <- 2
data$Y[which(data$Y == 'G2/M')] <- 3

X <- data[,2:(dim(data)[2]-1)]
X_labels <- data[,dim(data)[2]]
rownames(X)<-data[,1]
rownames(X)<-gsub("\\.","-",rownames(X))

# # SUBSET OF ENTRIES, JUST FOR TESTING
# X <- X[430:450,1:5]
# X_labels <- X_labels[430:450]

#Variance filtering
# varFiltering <- TRUE
# nvars <- 100
# if(varFiltering){
#   vars <- apply(X,2,var)
#   topfeatures <- names(sort(vars, dec=TRUE))[1:nvars]
#   X <- X[,topfeatures]
# }

#### Choose the starting and the ending points
boundary_init <- spathialBoundaryIds(X, X_labels, mode=1, from=2, to=1)
boundary_ids <- boundary_init$boundary_ids
X <- boundary_init$X
X_labels <- boundary_init$X_labels

#### Compute spathial (without filtering)
spathial_res <- spathialWay(X, X_labels, boundary_ids, NC, neighbors = 1)

# Analize the results
#Labels for each waypoint with knn
ppath_labels <- spathialLabels(X, X_labels, spathial_res)
#Plot the path in 2D using Rtsne
spathialPlot2D(X, X_labels, boundary_ids, spathial_res, perplexity_value=30)


# #### Compute spathial (with filtering)
# # Prefilter data
# filter_res <- spathialPrefiltering(X, X_labels, boundary_ids)
#
# X_filtered <- filter_res$X_filtered
# X_labels_filtered <- filter_res$X_labels_filtered
# X_garbage <- filter_res$X_garbage
# X_labels_garbage <- filter_res$X_labels_garbage
# boundary_ids_filtered <- filter_res$boundary_ids
#
# # Compute spathial
# spathial_res <- spathialWay(X_filtered, X_labels_filtered, boundary_ids_filtered, NC, neighbors = 1)
# message(length(spathial_res))
# message(spathial_res$perturbed_path)
# # Analize the results
# #Labels for each waypoint with knn
# ppath_labels <- spathialLabels(X_filtered, X_labels_filtered, spathial_res)
# #Plot the path in 2D using Rtsne
# spathialPlot2D(X_filtered, X_labels_filtered, boundary_ids, spathial_res, perplexity_value=30, X_garbage, X_labels_garbage)

####Correlation along the path
corr_values_principal_path <- spathialCorrelation(spathial_res)
