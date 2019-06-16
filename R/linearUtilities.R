## Initialize NC medoids with init_type rational.
#
#Args:
#[ndarray float] X: data matrix
#[int] n: number of medoids to be selected
#[string] init_type: rational to be used ('uniform' = randomly selected with uniform distribution, 'kpp' = k-means++ algorithm)
#[ndarray int] exclude_ids: blacklisted ids that shouldn't be selected
#
#Returns:
#[ndarray int] med_ids: indices of the medoids selected
#
initMedoids <- function(X, n, init_type, exclude_ids){
  N<-nrow(X)
  D<-ncol(X)
  med_ids<-rep(-1,n)
  if(init_type=="kpp"){
    # Select a random point (not a boundary id)
    med_ids[1]<-sample(setdiff(1:N,exclude_ids),1)

    # Now fill the path with n randomly picked medoids
    boundary_coords<-X[boundary_ids,]
    for(i in 1:(n-1)){
      # Calculate all distances between points up to here and boundaries
      Xmed_dst<-pdist(
        X,
        rbind(X[med_ids[1:i],],boundary_coords)
      )
      Xmed_dst<-as.matrix(Xmed_dst)^2
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
      med_id<-which.max(D2*D2_n)

      med_ids[i+1]<-med_id
    }
  } else if (init_type=='uniform'){
    med_ids<-sample(1:N,n)
  } else {
    stop("init_type not recognized.")
  }
  return(med_ids)
}

## Find the elbow in a fuction f, as the point on f with max distance
## from the line connecting f[0,:] and f[-1,:]
#
# Args:
# [ndarray float] f: function(Nx2 array in the form [x, f(x)])
#
# Returns:
# [int] elb_id: index of the elbow
#
find_elbow <- function(f){

}
