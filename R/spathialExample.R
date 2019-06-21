# Set Random Seed
set.seed(1)
# Number of optimization variables i.e. path waypoints
NC<-50

myfile<-system.file("extdata", "music_dataset.csv", package = "spathial")
data<-read.csv(myfile,as.is=TRUE,header=FALSE)
X <- data[,3:dim(data)[2]]
X_labels <- data[,1]

#X <- X[which(X_labels == 5 | X_labels == 42 | X_labels == 48), ]
#X_labels <- X_labels[which(X_labels == 5 | X_labels == 42 | X_labels == 48)]

X <- X[which(X_labels == 2 | X_labels == 48 | X_labels == 15), ]
X_labels <- X_labels[which(X_labels == 2 | X_labels == 48 | X_labels == 15)]

#X <- X[1:10,c("V10", "V14", "V17", "V22")]
#X_labels <- X_labels[1:10]

# Set row names
rownames(X)<-paste0("sam",rownames(X))
boundary_init <- spathial_boundary_ids(X, X_labels, mode=0, from=2, to=15)
boundary_ids <- boundary_init$boundary_ids
X <- boundary_init$X
X_labels <- boundary_init$X_labels

ppath <- spathialWay(X, boundary_ids, NC, prefiltering=FALSE)

spathial_2D_plot(X, X_labels, ppath)
