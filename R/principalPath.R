#' Minimizer
#'
#' Regularized K-mens for principal path, MINIMIZER.
#'
#' @param [ndarray float] X: data matrix
#' @param [ndarray float] init_W: initial waypoints matrix
#' @param [float] s: regularization parameter
#' @return [ndarray float] W: final waypoints matrix
#' @return [ndarray int] labels: final
#' @
rkm <- function(X, init_W, s){

}

#' Prefiltering
#'
#' Regularized K-mens for principal path, PREFILTER.
#'
#' @param [ndarray float] X: data matrix
#' @param [ndarray int] boundary_ids: start/end waypoints as sample indices
#' @param [int] Nf: number of filter centroids
#' @param [int] k: number of nearest neighbor for the penalized graph
#' @param [float] p: penalty factor for the penalized graph
#' @param [float] T: filter threshold
#' @return [ndarray float] X_filtered
#' @return [ndarray int] X_labels_filtered
#' @return [ndarray int] boundary_ids_filtered
#' @return [ndarray float] X_garbage
#' @return [ndarray int] X_labels_garbage
rkm_prefilter <- function(X, boundary_ids, Nf=200, k=5, p=1000, T=0.1){

}

#' Model Selection
#'
#' Regularized K-menans for principal path, MODEL SELECTION, variance on waypoints interdistance.
#'
#' @param [ndarray float] models: matrix with path models, shape, N_models x N x (NC+2)
#' @param [ndarray float] s_span: array with values of the reg parameter for each model (sorted in decreasing order, with 0 as last value)
#' @param [ndarray float] X: data matrix
#' @return [ndarray float] W_dst_var: array with values of variance for each model
rkm_MS_pathvar <- function(models, s_span, X){

}
