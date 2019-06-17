# ### NOTE: this is a stub, the example should go in the vignette at the end
# library(spathial)

### Initialize rkm_2D example
# Set Random Seed
set.seed(1)

# Number of optimization variables i.e. path waypoints
NC<-50

# Read 2D input file (constellation)
myfile<-system.file("extdata", "2D_constellation.csv", package = "spathial")
X<-read.csv(myfile,as.is=TRUE,header=FALSE)

# NOTE: A safety measure to prevent using integers as character indexes
rownames(X)<-paste0("sam",rownames(X))

### Select Boundaries
plot(X,pch=20,col="black",main="Click to select path start and end points")
if(FALSE){ # Set to TRUE for manual selection
  boundary_ids<-rownames(X)[identify(X,n=2,plot=FALSE)]
}
boundary_ids<-rownames(X)[c(424,593)]
points(
  X[boundary_ids,],pch="x",col="red",cex=4,
  xlab="Dimension 1",ylab="Dimension 2"
)

### Prefilter the data (function pp.rkm_prefilter)
prefiltered<-rkm_prefilter(X,boundary_ids,plot_ax=TRUE)
X<-prefiltered$X_filtered
boundary_ids<-prefiltered$boundary_ids_filtered
X_g<-prefiltered$X_garbage
rm(prefiltered)


### Initialize waypoints
waypoint_ids<-initMedoids(X, NC, 'kpp', boundary_ids)
waypoint_ids<-c(boundary_ids[1],waypoint_ids,boundary_ids[2])
init_W<-X[waypoint_ids,]


### Annealing with rkm
s_span<-pracma::logspace(5,-5,n=NC)
s_span<-c(s_span,0)
#models<-array(data=NA,dim=c(length(s_span),NC+2,ncol(X)))
#s<-s_span[1]
models<-list()

for(i in 1:length(s_span)){
  s<-s_span[i]
  W<-rkm(X,init_W,s,plot_ax=TRUE)
  init_W<-W
  models[[as.character(s)]]<-W
  #models[i,,]<-W
}


### Interactive model selection




