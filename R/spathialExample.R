### NOTE: this is a stub, the example should go in the vignett at the end

### Initialize rkm_2D example
# Set Random Seed
set.seed(1)
# Number of optimization variables i.e. path waypoints
nc<-50
# Read 2D input file (constellation)
myfile<-system.file("extdata", "2D_constellation.csv", package = "spathial")
X<-read.csv(myfile,as.is=TRUE,header=FALSE)
N=nrow(X)
d=ncol(X)

### Select Boundaries
# NOTE: the example could benefit from a manually selected pair
# boundary_ids=lu.getMouseSamples2D(X,2)
boundary_ids=c(424,593)
plot(X,pch=20,col="black")
points(
  X[boundary_ids,],pch="x",col="red",cex=4,
  xlab="Dimension 1",ylab="Dimension 2"
)

### Prefilter the data (function pp.rkm_prefilter)
Nf<-200
med_ids<-initMedoids(X,n=Nf-2,init_type="kpp",exclude_ids=boundary_ids)
med_ids<-c(boundary_ids[1],med_ids,boundary_ids[2])






