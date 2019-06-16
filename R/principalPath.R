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
#
# Returns:
#
# [ndarray float] X_filtered
# [ndarray int] X_labels_filtered
# [ndarray int] boundary_ids_filtered
# [ndarray float] X_garbage
# [ndarray int] X_labels_garbage
rkm_prefilter <- function(X, boundary_ids, Nf=200, k=5, p=1000, T=0.1){
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
  medmed_dst_p[1,Nf]<-0
  medmed_dst_p[Nf,1]<-0
  # Find shortest path using dijkstra



}

#Model Selection
#Regularized K-menans for principal path, MODEL SELECTION, variance on waypoints interdistance.
#[ndarray float] models: matrix with path models, shape, N_models x N x (NC+2)
#[ndarray float] s_span: array with values of the reg parameter for each model (sorted in decreasing order, with 0 as last value)
#[ndarray float] X: data matrix
#[ndarray float] W_dst_var: array with values of variance for each model
rkm_MS_pathvar <- function(models, s_span, X){

}
