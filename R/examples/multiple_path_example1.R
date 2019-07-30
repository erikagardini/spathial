<<<<<<< HEAD:R/examples/multiple_path_example1.R
####Load input data
myfile<-system.file("extdata", "liver_tcga_example1.csv", package = "spathial")
data<-read.csv(myfile,as.is=TRUE,header=TRUE)
X <- data[,2:(ncol(data)-1)]
rownames(X)<-paste0("sam",rownames(X))
X_labels <- data$Y

#### Choose the starting and the ending points
boundary_init <- spathialBoundaryIds(X, X_labels, mode=2, from=2, to=1)
boundary_ids <- boundary_init$boundary_ids
X <- boundary_init$X
X_labels <- boundary_init$X_labels

# Compute spathial
NC <- 50
spathial_res <- spathialWay(X, X_labels, boundary_ids, NC, neighbors = 1)

#Plot the path in 2D using Rtsne
spathialPlot(X, X_labels, boundary_ids, spathial_res, perplexity_value=30)

# ####Correlation along the path
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
=======
# ####Load input data
# myfile<-system.file("extdata", "liver_tcga_example1.csv", package = "spathial")
# data<-read.csv(myfile,as.is=TRUE,header=TRUE)
# X <- data[,2:(ncol(data)-1)]
# rownames(X)<-paste0("sam",rownames(X))
# X_labels <- data$Y
#
# X <- X[,1:5]
#
# #### Choose the starting and the ending points
# boundary_init <- spathialBoundaryIds(X, X_labels, mode=2, from=2, to=1)
# boundary_ids <- boundary_init$boundary_ids
# X <- boundary_init$X
# X_labels <- boundary_init$X_labels
#
# # Compute spathial
# NC <- 50
# spathial_res <- spathialWay(X, X_labels, boundary_ids, NC, neighbors = 2)
#
# #Plot the path in 2D using Rtsne
# spathialPlot(X, X_labels, boundary_ids, spathial_res, perplexity_value=30)
#
# # ####Correlation along the path
# statistics <- spathialStatistics(spathial_res)
#
# #Plot the results
# corr_values <- as.numeric(unlist(statistics$correlations))
#
# sub_correlation_matrix <- matrix(data=0, nrow = (length(corr_values)), ncol = 1)
# colnames(sub_correlation_matrix) <- "pp"
# rownames(sub_correlation_matrix) <- colnames(X)
# sub_correlation_matrix[,ncol(sub_correlation_matrix)] <- corr_values
#
# library(corrplot)
# library(RColorBrewer)
# corrplot(sub_correlation_matrix, method="circle", is.corr=FALSE, cl.pos = "r", cl.ratio = 1, cl.lim = c(-1,1))
>>>>>>> 85d5d26208eb45d425fb2fd54e009c33471852e2:R/multiple_path_example.R
