# Set Random Seed
set.seed(1)
# Number of optimization variables i.e. path waypoints
NC <- 2

myfile<-system.file("extdata", "thyroid_tcga.csv", package = "spathial")
data<-read.csv(myfile,as.is=TRUE,header=TRUE)
X <- data[,2:(dim(data)[2]-1)]
X_labels <- data[,dim(data)[2]]
# Set row names
rownames(X)<-data[,1]

X <- X[,1:5]
#Variance filtering
#varFiltering <- TRUE
#nvars <- 5
#if(varFiltering){
 # vars <- apply(X,2,var)
  #topfeatures <- names(sort(vars, dec=TRUE))[1:nvars]
 # X <- X[,topfeatures]
#}

X <- X[438:445,]
X_labels <- X_labels[438:445]
#Choose the starting and the ending points
message("BOUNDARY INIT")
boundary_init <- spathial_boundary_ids(X, X_labels, mode=1, from=2, to=1)
boundary_ids <- boundary_init$boundary_ids
X <- boundary_init$X
X_labels <- boundary_init$X_labels
message(boundary_ids)
message(dim(X))
message(dim(X_labels))

#Compute spathial
message("COMPUTE SPATHIAL")
spathial_res <- spathial_way_multiple(X, X_labels, boundary_ids, NC, prefiltering=FALSE, negb = 2)

#Labels for each waypoint with knn
message("SPATHIAL LABELS")
ppath_labels <- spathial_labels(X, X_labels, spathial_res$ppath)

#Plot the path in 2D using tsne
message("SPATHIAL PLOT")
spathial_2D_plot(X, X_labels, boundary_ids, spathial_res$ppath)

#Correlation along the path
#message("SPATHIAL CORRELATION")
#corr_values <- spathial_corr(spathial_res)
#message(dim(corr_values))
