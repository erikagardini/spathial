#Regularized K-mens for principal path, MINIMIZER.
#[ndarray float] X: data matrix
#[ndarray float] init_W: initial waypoints matrix
#[float] s: regularization parameter
#[ndarray float] W: final waypoints matrix
#[ndarray int] labels: final
rkm <- function(X, init_W, s){

}

## Regularized K-means for principal path, PREFILTER.
#
# Args:
# [ndarray float] X: data matrix
# [ndarray int] boundary_ids: start/end waypoints as sample indices
# [int] Nf: number of filter centroids
# [int] k: number of nearest neighbor for the penalized graph
# [float] p: penalty factor for the penalized graph
# [float] T: filter threshold
# plot_ax: boolean, whether the post-filtering graph should be plotted
#
# Returns:
#
# [ndarray float] X_filtered
# [ndarray int] X_labels_filtered
# [ndarray int] boundary_ids_filtered
# [ndarray float] X_garbage
# [ndarray int] X_labels_garbage

#' rkm_prefilter
#'
#' Regularized K-means for principal path: prefiltering
#'
#' @param X a data matrix
#' @param boundary_ids
#' @param Nf
#' @param k
#' @param T
#' @param plot_ax boolean, whether the post-filtering graph should be plotted
#' @return A list with letters and numbers.
#' \itemize{
#'   \item X_filtered - The filtered data matrix
#'   \item boundary_ids_filtered - The filtered boundary ids
#'   \item X_garbage - The data matrix that was discarded
#' }
#' @export
rkm_prefilter <- function(X, boundary_ids, Nf=200, k=5, p=1000, T=0.1,
                          plot_ax=FALSE){
  N=nrow(X)
  d=ncol(X)
  # Pick Nf medoids with k-means++ and compute pairwise distance matrix
  med_ids<-initMedoids(X,n=Nf-2,init_type="kpp",exclude_ids=boundary_ids)
  med_ids<-c(boundary_ids[1],med_ids,boundary_ids[2])
  medmed_dst<-as.matrix(dist(X[med_ids,],method="euclidean"))^2

  # Build k-nearest-neighbor penalized matrix
  knn_ids<-t(apply(medmed_dst,1,order))
  medmed_dst_p<-medmed_dst*p
  for (i in 1:Nf){
    for(j in 1:k){
      k<-knn_ids[i,j]
      medmed_dst_p[i,k]<-medmed_dst[i,k]
      medmed_dst_p[k,i]<-medmed_dst[k,i]
    }
  }
  medmed_dst_p[1,Nf]<-0 # Not sure why this is done
  medmed_dst_p[Nf,1]<-0 # Not sure why this is done

  # Find shortest path using dijkstra
  # NOTE: I didn't calculate the predecessors and simply got the shortest path
  g<-graph.adjacency(medmed_dst_p, weighted=TRUE)
  # path_dst<-igraph::distances(g,v=1,algorithm="dijkstra")
  path<-igraph::shortest_paths(g,from=1,to=ncol(medmed_dst_p))
  path<-names(path$vpath[[1]])

  # Filter out medoids too close to the shortest path
  T<-T*mean(medmed_dst)
  to_filter_ids<-c()
  for(i in path){
    to_filter_ids<-c(names(which(medmed_dst[i,]<T,useNames=TRUE)),to_filter_ids)
  }
  to_filter_ids<-unique(setdiff(to_filter_ids,path))
  to_keep_ids<-setdiff(colnames(medmed_dst),to_filter_ids)
  Xmed_dst<-pdist::pdist(
    X,
    X[to_keep_ids,]
  )
  Xmed_dst<-as.matrix(Xmed_dst)^2
  rownames(Xmed_dst)<-rownames(X)
  colnames(Xmed_dst)<-to_keep_ids
  u<-names(apply(Xmed_dst,1,which.min))
  filter_mask<-setNames(rep(FALSE,N),rownames(X))
  for(i in 1:N){
    if(u[i]%in%path){
      filter_mask[i]<-TRUE
    }
  }

  # Convert boundary indices
  boundary_ids_filtered<-boundary_ids
  boundary_ids_filtered[1]<-boundary_ids[1] - boundary_ids[1] + sum(filter_mask[1:boundary_ids[1]])
  boundary_ids_filtered[2]<-boundary_ids[2] - boundary_ids[2] + sum(filter_mask[1:boundary_ids[2]])

  # Plot filter figure
  if(plot_ax){
    xrange<-abs(max(X[!filter_mask,1])-min(X[!filter_mask,1]))
    plot(X[!filter_mask,1],X[!filter_mask,2],
         xlab="Dimension 1",ylab="Dimension 2",
         main="Post-filtering Path plot",col="orange",pch=1,
         xlim=c(min(X[!filter_mask,1]),max(X[!filter_mask,1])+xrange/2)) # data filtered out
    points(X[filter_mask,1],X[filter_mask,2],pch=19,col="blue") # data kept NOTE: THIS ONE IS WRONG, recheck filter_mask
    points(X[med_ids,1],X[med_ids,2],col="red",pch=1) # filter medoids
    points(X[to_filter_ids, 1], X[to_filter_ids, 2],pch=4,cex=1.2) # filter medoids dropped
    lines(X[path,1], X[path,2],col="darkgreen",lwd=3,pch=19,type="o") # filter shortest path
    text(X[filter_mask,][boundary_ids_filtered,1],
         X[filter_mask,][boundary_ids_filtered,2],
         labels=c("B1","B2"),font=2,cex=1.2) # Boundary samples
    legend("right",
           legend=c(
             "data filtered out",
             "data kept",
             "filter medoids",
             "filter medoids dropped",
             "filter shortest path",
             "boundary samples"
           ),bg="#FFFFFF22",
           pch=c(1,19,1,4,19,66),
           col=c("orange","blue","red","black","darkgreen","black"),
           lty=c(0,0,0,0,1,0),lwd=c(0,0,0,0,2,0)
    )
  }
  # Return ouput
  outlist<-list(X[filter_mask,], boundary_ids_filtered, X[!filter_mask,])
  return(outlist)
}

#Model Selection
#Regularized K-menans for principal path, MODEL SELECTION, variance on waypoints interdistance.
#[ndarray float] models: matrix with path models, shape, N_models x N x (NC+2)
#[ndarray float] s_span: array with values of the reg parameter for each model (sorted in decreasing order, with 0 as last value)
#[ndarray float] X: data matrix
#[ndarray float] W_dst_var: array with values of variance for each model
rkm_MS_pathvar <- function(models, s_span, X){

}
