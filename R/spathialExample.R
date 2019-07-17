####Load input data
myfile<-system.file("extdata", "2D_constellation.csv", package = "spathial")
data<-read.csv(myfile,as.is=TRUE,header=TRUE)
X <- data[,2:(ncol(data)-1)]
rownames(X)<-paste0("sam",rownames(X))
X_labels <- data$Y

colors <- rainbow(length(table(X_labels)))
colors_labels <- sapply(X_labels, function(x){colors[x]})
plot(X[,1],X[,2], col=colors_labels, pch=as.character(X_labels))

#### Choose the starting and the ending points
boundary_init <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
boundary_ids <- boundary_init$boundary_ids
X <- boundary_init$X
X_labels <- boundary_init$X_labels

#Plot the results
boundaries <- X[which(rownames(X) == boundary_ids[1] | rownames(X) == boundary_ids[2]),]
plot(X[,1],X[,2], col=colors_labels, pch=as.character(X_labels))
points(boundaries[,1],boundaries[,2], pch="x",col="black",cex=4)

#### Compute spathial (with filtering)
# Prefilter data
filter_res <- spathialPrefiltering(X, X_labels, boundary_ids)
mask <- filter_res$mask
boundary_ids <- filter_res$boundary_ids

#Plot the results
boundaries <- X[which(rownames(X) == boundary_ids[1] | rownames(X) == boundary_ids[2]),]
X_garbage <- X[!mask,]

plot(X[,1],X[,2], col=colors_labels, pch=as.character(X_labels))
points(boundaries[,1],boundaries[,2], pch="x",col="black",cex=4)
points(X_garbage[,1],X_garbage[,2], col="gray", pch=4)

# Compute spathial
X_filtered <- X[mask,]
X_labels_filtered <- X_labels[mask]
NC <- 50
spathial_res_filtered <- spathialWay(X_filtered, X_labels_filtered, boundary_ids, NC, neighbors = 1)
spathial_res <- spathialWay(X, X_labels, boundary_ids, NC, neighbors = 1)

#Plot the path in 2D using Rtsne
spathialPlot(X, X_labels, boundary_ids, spathial_res_filtered, perplexity_value=30, mask)
spathialPlot(X, X_labels, boundary_ids, spathial_res, perplexity_value=30)

# Analize the results
#Labels for each waypoint with knn
ppath_labels_filtered <- spathialLabels(X_filtered, X_labels_filtered, spathial_res_filtered)
ppath_labels <- spathialLabels(X, X_labels, spathial_res)

#Plot the results
ppath_labels_filtered <- as.vector(ppath_labels_filtered)
colors_labels_ppath <- sapply(ppath_labels_filtered, function(y){colors[as.integer(y)]})
plot(c(1:length(ppath_labels_filtered)), c(ppath_labels_filtered), col=colors_labels_ppath, pch=as.character(ppath_labels_filtered))

ppath_labels <- as.vector(ppath_labels)
colors_labels_ppath <- sapply(ppath_labels, function(y){colors[as.integer(y)]})
plot(c(1:length(ppath_labels)), c(ppath_labels), col=colors_labels_ppath, pch=as.character(ppath_labels))

# ####Correlation along the path
corr_values_principal_path <- spathialStatistics(spathial_res)

#Plot the results
corr_values <- as.numeric(unlist(corr_values_principal_path$correlations))
correlation_matrix <- matrix(data=0, nrow = (length(corr_values)+1), ncol = (length(corr_values)+1))
colnames(correlation_matrix) <- c(colnames(X), "pp")
rownames(correlation_matrix) <- c(colnames(X), "pp")
correlation_matrix[,ncol(correlation_matrix)] <- c(corr_values, 0)
correlation_matrix[nrow(correlation_matrix),] <- c(corr_values, 0)

library(corrplot)
library(RColorBrewer)
corrplot(correlation_matrix, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))
