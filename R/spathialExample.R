# Set Random Seed
set.seed(1)
# Number of optimization variables i.e. path waypoints
NC <- 50


####Init input data
myfile<-system.file("extdata", "example.csv", package = "spathial")
data<-read.csv(myfile,as.is=TRUE,header=TRUE)
X <- data[,2:(dim(data)[2]-1)]
rownames(X)<-paste0("sam",rownames(X))
X_labels <- data$Y

#Plot the results
set.seed(1)
tsne_res <- Rtsne::Rtsne(as.matrix(X), dims = 2, perplexity = 30)
points_2D <- tsne_res$Y
boundary_ids_2D <- points_2D[which(rownames(X) == boundary_ids[1] | rownames(X) == boundary_ids[2]),]

plot(points_2D[,1],points_2D[,2], col=X_labels, pch=19)

#### Choose the starting and the ending points
boundary_init <- spathialBoundaryIds(X, X_labels, mode=2, from=1, to=2)
boundary_ids <- boundary_init$boundary_ids
X <- boundary_init$X
X_labels <- boundary_init$X_labels

#Plot the results
set.seed(1)
tsne_res <- Rtsne::Rtsne(as.matrix(X), dims = 2, perplexity = 30)
points_2D <- tsne_res$Y
boundary_ids_2D <- points_2D[which(rownames(X) == boundary_ids[1] | rownames(X) == boundary_ids[2]),]

plot(points_2D[,1],points_2D[,2], col=X_labels, pch=19)
points(boundary_ids_2D[,1],boundary_ids_2D[,2], pch="x",col="green",cex=4)

#### Compute spathial (without filtering)
spathial_res <- spathialWay(X, X_labels, boundary_ids, NC, neighbors = 1)

# Analize the results
#Labels for each waypoint with knn
ppath_labels <- spathialLabels(X, X_labels, spathial_res)
#Plot the path in 2D using Rtsne
spathialPlot2D(X, X_labels, boundary_ids, spathial_res, perplexity_value=30)


#### Compute spathial (with filtering)
# Prefilter data
filter_res <- spathialPrefiltering(X, X_labels, boundary_ids)

mask <- filter_res$mask
boundary_ids <- filter_res$boundary_ids

#Plot the results
set.seed(1)
tsne_res <- Rtsne::Rtsne(as.matrix(X), dims = 2, perplexity = 30)
points_2D <- tsne_res$Y

boundary_ids_2D <- points_2D[which(rownames(X) == boundary_ids[1] | rownames(X) == boundary_ids[2]),]
X_garbage_2D <- points_2D[!mask,]

plot(points_2D[,1],points_2D[,2], col=X_labels, pch=19)
points(boundary_ids_2D[,1],boundary_ids_2D[,2], pch="x",col="green",cex=4)
points(X_garbage_2D[,1],X_garbage_2D[,2], col="gray", pch=4)



# Compute spathial
spathial_res <- spathialWay(X_filtered, X_labels_filtered, boundary_ids_filtered, NC, neighbors = 1)
message(length(spathial_res))
message(spathial_res$perturbed_path)
# Analize the results
#Labels for each waypoint with knn
ppath_labels <- spathialLabels(X_filtered, X_labels_filtered, spathial_res)
#Plot the path in 2D using Rtsne
spathialPlot2D(X_filtered, X_labels_filtered, boundary_ids, spathial_res, perplexity_value=30, X_garbage, X_labels_garbage)
#
# ####Correlation along the path
# corr_values_principal_path <- spathialCorrelation(spathial_res)
