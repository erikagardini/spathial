####Load input data
myfile<-system.file("extdata", "2D_constellation.csv", package = "spathial")
data<-read.csv(myfile,as.is=TRUE,header=TRUE)
X <- data[,2:(ncol(data)-1)]
rownames(X)<-paste0("sam",rownames(X))
X_labels <- data$Y
plot(X[,1],X[,2], col=X_labels, pch=X_labels)

#### Choose the starting and the ending points
boundary_init <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
boundary_ids <- boundary_init$boundary_ids
X <- boundary_init$X
X_labels <- boundary_init$X_labels

#Plot the results
boundaries <- X[which(rownames(X) == boundary_ids[1] | rownames(X) == boundary_ids[2]),]
plot(X[,1],X[,2], col=X_labels, pch=X_labels)
points(boundaries[,1],boundaries[,2], pch="x",col="red",cex=4)

#### Compute spathial (with filtering)
# Prefilter data
filter_res <- spathialPrefiltering(X, X_labels, boundary_ids)
mask <- filter_res$mask
boundary_ids <- filter_res$boundary_ids

#Plot the results
boundary_ids_2D <- X[which(rownames(X) == boundary_ids[1] | rownames(X) == boundary_ids[2]),]
X_garbage_2D <- X[!mask,]

plot(X[,1],X[,2], col=X_labels, pch=X_labels)
points(boundary_ids_2D[,1],boundary_ids_2D[,2], pch="x",col="red",cex=4)
points(X_garbage_2D[,1],X_garbage_2D[,2], col="blue", pch=4)

# Compute spathial
X_filtered <- X[mask,]
X_labels_filtered <- X_labels[mask]
NC <- 50
spathial_res <- spathialWay(X, X_labels, boundary_ids, NC, neighbors = 1)

# Analize the results
#Labels for each waypoint with knn
ppath_labels <- spathialLabels(X, X_labels, spathial_res)
#Plot the path in 2D using Rtsne
spathialPlot2D(X, X_labels, boundary_ids, spathial_res, perplexity_value=30)#, mask)
#
# ####Correlation along the path
# corr_values_principal_path <- spathialCorrelation(spathial_res)
