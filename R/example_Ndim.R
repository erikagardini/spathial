#Load data
myfile<-system.file("extdata", "thyroid_tcga_example1.csv", package = "spathial")
data<-read.csv(myfile,as.is=TRUE,header=TRUE)
X <- data[,2:(dim(data)[2]-1)]
X_labels <- data[,dim(data)[2]]
rownames(X)<-data[,1]

# Choose the starting and the ending points
boundary_init <- spathialBoundaryIds(X, X_labels, mode=2, from=2, to=1)#3, from="TCGA-DD-A39W-11A-11R-A213-07", to="TCGA-G3-AAV2-01A-11R-A37K-07")
boundary_ids <- boundary_init$boundary_ids
X <- boundary_init$X
X_labels <- boundary_init$X_labels

# Compute spathial
NC <- 50
spathial_res <- spathialWay(X, X_labels, boundary_ids, NC, neighbors = 1)

# Analize the results
#Labels for each waypoint with knn
ppath_labels <- spathialLabels(X, X_labels, spathial_res)
ppath_labels <- as.vector(ppath_labels)
colors_labels_ppath <- sapply(ppath_labels, function(y){colors[as.integer(y)]})
plot(c(1:length(ppath_labels)), c(ppath_labels), col=colors_labels_ppath, pch=as.character(ppath_labels))

#Plot the path in 2D using Rtsne
spathialPlot(X, X_labels, boundary_ids, spathial_res, perplexity_value=30)

#Correlation along the path
statistics <- spathialStatistics(spathial_res)

#Plot the results
corr_values <- as.numeric(unlist(statistics$correlations))
arg_sort <- rank(-abs(corr_values))
indexes <- which(arg_sort < 20)

sub_correlation_matrix <- matrix(data=0, nrow = (length(indexes)), ncol = 1)
colnames(sub_correlation_matrix) <- "pp"
rownames(sub_correlation_matrix) <- colnames(X[indexes])
sub_correlation_matrix[,ncol(sub_correlation_matrix)] <- corr_values[indexes]

library(corrplot)
library(RColorBrewer)
corrplot(sub_correlation_matrix, method="circle", is.corr=FALSE, cl.pos = "r", cl.ratio = 1, cl.lim = c(-1,1))
