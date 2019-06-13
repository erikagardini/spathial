#' Compute principal path
#'
#' Get the coordinates of the waypoints of the principal path
#'
#' @param [ndarray float] X: data points
#' @param [ndarray int] boundaries: starting and ending points
#' @param [int] NC: number of waypoints
#' @param [ndarray logic] prefiltering: TRUE/FALSE
#' @return [ndarray float] ppath: waypoints
#' @export
spathial <- function(X, boundaries, NC, prefiltering){

}

#' Find labels
#'
#' Get the label of each waypoint accordin to the neighbourhood
#'
#'@param [ndarray float] X: data points
#'@param [ndarray int] X_labels: labels of the data points
#'@param [ndarray float] ppath: waypoints
#'@return [ndarray int] ppath_labels: labels of the waypoints
#'@export
spathial_labels <- function(X, X_labels, ppath){

}

#' 2D spathial
#'
#' Get the 2D coordinates of each waypoint (using t-SNE algorithm for the dimensionality reduction)
#'
#' @param [ndarray float] ppath: waypoints
#' @return [ndarray float] 2D_ppath: 2D coordinates of the waypoints
#' @export
spathial_2D <- function(ppath){

}

#CORRELATION?????
