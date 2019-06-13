#' Initialize medoids
#'
#' Initialize NC medoids with init_type rational.
#'
#' @param [ndarray float] X: data matrix
#' @param [int] n: number of medoids to be selected
#' @param [string] init_type: rational to be used ('uniform' = randomly selected with uniform distribution, 'kpp' = k-means++ algorithm)
#' @param [ndarray int] exclude_ids: blacklisted ids that shouldn't be selected
#' @return [ndarray int] med_ids: indices of the medoids selected
initMedoids <- function(X, n, init_type, exclude_ids){

}

#' Find the elbow id
#'
#' Find the elbow in a fuction f, as the point on f with max distance from the line connecting f[0,:] and f[-1,:]
#'
#' @param [ndarray float] f: function(Nx2 array in the form [x, f(x)])
#' @return [int] elb_id: index of the elbow
find_elbow <- function(f){

}

prova <- function(){
  x <- 2
  return(x)
}

