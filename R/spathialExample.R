# ### NOTE: this is a stub, the example should go in the vignette at the end
# library(spathial)

### Initialize rkm_2D example
# Set Random Seed
set.seed(1)

# Number of optimization variables i.e. path waypoints
NC<-50

# Read 2D input file (constellation)
#myfile<-system.file("extdata", "2D_constellation.csv", package = "spathial")
#X<-read.csv(file="/Users/erikagardini/spathial/inst/extdata/music_dataset.csv",as.is=TRUE,header=FALSE)
myfile<-system.file("extdata", "music_dataset.csv", package = "spathial")
X<-read.csv(myfile,as.is=TRUE,header=FALSE)
X <- X[,3:dim(X)[2]]
# NOTE: A safety measure to prevent using integers as character indexes
rownames(X)<-paste0("sam",rownames(X))

### Select Boundaries
plot(X[,57],X[,501], pch=20,col="black",main="Click to select path start and end points")
if(FALSE){ # Set to TRUE for manual selection
  boundary_ids<-rownames(X)[identify(X,n=2,plot=FALSE)]
}
boundary_ids<-rownames(X)[c(21,579)]
points(
  X[boundary_ids,57], X[boundary_ids,501],pch="x",col="red",cex=4,
  xlab="Dimension 1",ylab="Dimension 2"
)

#waypoint_ids<-initMedoids(X, NC, 'kpp', boundary_ids)
pp <- spathialWay(X, boundary_ids, NC, FALSE)

pp_2D <- spathial_2D(pp)
plot(X[,57],X[,501])
lines(pp_2D$Y[,1], pp_2D$Y[,2],lwd=3,col="red",type="o",pch=15)
#id <- spathial(X, boundary_ids, NC, TRUE)
#plot(waypoints,pch=20,col="black",main="Waypoints")
