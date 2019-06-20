# Set Random Seed
set.seed(1)
# Number of optimization variables i.e. path waypoints
NC<-50

myfile<-system.file("extdata", "music_dataset.csv", package = "spathial")
data<-read.csv(myfile,as.is=TRUE,header=FALSE)
X <- data[,3:dim(data)[2]]
X_labels <- data[,1]

# Set row names
rownames(X)<-paste0("sam",rownames(X))
boundary_init <- spathial_boundary_ids(X, X_labels, mode=0, from=2, to=48)
boundary_ids <- boundary_init$boundary_ids
X <- boundary_init$X
X_labels <- boundary_init$X_labels


X_2D <- spathial_2D(X)

pp <- spathialWay(X, boundary_ids, NC, prefiltering=FALSE)


pp_2D <- spathial_2D(pp)
#plot(X_2D$Y[,1],X_2D$Y[,2])
#lines(pp_2D$Y[,1], pp_2D$Y[,2],lwd=3,col="red",type="o",pch=15)
