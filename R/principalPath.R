#Regularized K-mens for principal path, MINIMIZER.
#[ndarray float] X: data matrix
#[ndarray float] init_W: initial waypoints matrix
#[float] s: regularization parameter
#[ndarray float] W: final waypoints matrix
#[ndarray int] labels: final
rkm <- function(X, init_W, s){

}

#Regularized K-mens for principal path, PREFILTER.
#[ndarray float] X: data matrix
#[ndarray int] boundary_ids: start/end waypoints as sample indices
#[int] Nf: number of filter centroids
#[int] k: number of nearest neighbor for the penalized graph
#[float] p: penalty factor for the penalized graph
#[float] T: filter threshold
#[ndarray float] X_filtered
#[ndarray int] X_labels_filtered
#[ndarray int] boundary_ids_filtered
#[ndarray float] X_garbage
#[ndarray int] X_labels_garbage
rkm_prefilter <- function(X, boundary_ids, Nf=200, k=5, p=1000, T=0.1){

}

#Model Selection
#Regularized K-menans for principal path, MODEL SELECTION, variance on waypoints interdistance.
#[ndarray float] models: matrix with path models, shape, N_models x N x (NC+2)
#[ndarray float] s_span: array with values of the reg parameter for each model (sorted in decreasing order, with 0 as last value)
#[ndarray float] X: data matrix
#[ndarray float] W_dst_var: array with values of variance for each model
rkm_MS_pathvar <- function(models, s_span, X){

}
