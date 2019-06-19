
#' Find labels
#'
#' Get the label of each waypoint accordin to the neighbourhood
#'
#' @param [ndarray float] X: data points
#' @param [ndarray int] X_labels: labels of the data points
#' @param [ndarray float] ppath: waypoints
#' @return [ndarray int] ppath_labels: labels of the waypoints
#' @export
spathial_labels <- function(X, X_labels, ppath){
  lbl <- knn(X, ppath, X_labels, k=1)
  return(lbl)
}

#' 2D spathial
#'
#' Get the 2D coordinates of each waypoint (using t-SNE algorithm for the dimensionality reduction)
#'
#' @param [ndarray float] ppath: waypoints
#' @return [ndarray float] 2D_ppath: 2D coordinates of the waypoints
#' @export
spathial_2D <- function(ppath){
  library(Rtsne)
  set.seed(1)
  ttt<-Rtsne(as.matrix(ppath), dims = 2, perplexity = 5)
  return(ttt)
}


