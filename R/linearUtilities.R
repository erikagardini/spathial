#' initMedoids
#'
#' Initialize NC medoids with init_type rational.
#'
#' @param X a data matrix
#' @param n integer number of medoids to be selected
#' @param init_type character, the rational to be used
#' ('uniform' = randomly selected with uniform distribution, 'kpp' = k-means++ algorithm)
#' @param boundary_ids blacklisted ids that shouldn't be selected
#' @return A vevtor of medoids selected
#' @export
initMedoids <- function(X, n, init_type, boundary_ids){
  N<-nrow(X)
  d<-ncol(X)
  med_ids<-rep(-1,n)
  if(init_type=="kpp"){
    # Select a random point (not a boundary id)
    med_ids[1]<-sample(setdiff(rownames(X),boundary_ids),1)

    # Now fill the path with n randomly picked medoids
    boundary_coords<-X[boundary_ids,]
    for(i in 1:(n-1)){
      # Calculate all distances between points up to here and boundaries
      Xmed_dst<-pracma::distmat(
        as.matrix(X),
        as.matrix(rbind(X[med_ids[1:i],],boundary_coords))
      )
      Xmed_dst<-Xmed_dst^2
      D2<-apply(Xmed_dst,1,min)
      D2_n<-1.0/sum(D2)

      ## NOTE: the following code is the original one but it contains some slow
      ## procedures in R. I rewrote something similar below. The idea is to
      ## select a distance with a high chance to be long, and to avoid taking
      ## the same medoid twice (there is likely an easier way to do this)
      # accepted<-FALSE
      # while(!accepted){
      #   med_id<-sample(1:N,1)
      #   rnumber<-runif(1)
      #   if(rnumber<(D2[med_id]*D2_n)){
      #     accepted<-TRUE
      #   }
      # }
      ## The following line summarizes the previous code without the
      ## randomness. If randomness is necessary, we will change this
      med_id<-rownames(X)[which.max(D2*D2_n)]

      med_ids[i+1]<-med_id
    }
  } else if (init_type=='uniform'){
    med_ids<-sample(rownames(X),n)
  } else {
    stop("init_type not recognized.")
  }
  return(med_ids)
}

#' Find the elbow in a fuction f, as the point on f with max distance
#' from the line connecting f[0,:] and f[-1,:]
#'
#' @param f: function(Nx2 array in the form [x, f(x)])
#' @return elb_id: index of the elbow
#' @export
find_elbow <- function(f){
  ps <- array(c(f[1,1], f[1,2]), dim = c(2))
  pe <- array(c(f[nrow(f),1], f[nrow(f),2]), dim = c(2))
  p_line_dst = array(data = 0.0, (nrow(f)-2))
  for(i in (2:(nrow(f)-1))){
    p <- array(c(f[i,1],f[i,2]), dim= c(2))

    mtx_1 <- pe-ps
    mtx_2 <- ps-p
    mtx_3 <- pe-ps
    cp <- ((mtx_1[1]*mtx_2[2]) - (mtx_1[2]*mtx_2[1]))
    den <- norm(as.matrix(cp), "F")
    num <- norm(as.matrix(mtx_3), "F")

    p_line_dst[i-1] <- (den / num)
  }
  elb_id = (which.max(p_line_dst)) + 1
  return(elb_id)
}
