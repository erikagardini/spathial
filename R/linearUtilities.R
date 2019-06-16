#Initialize NC medoids with init_type rational.
#Args:
#[ndarray float] X: data matrix
#[int] n: number of medoids to be selected
#[string] init_type: rational to be used ('uniform' = randomly selected with uniform distribution, 'kpp' = k-means++ algorithm)
#[ndarray int] exclude_ids: blacklisted ids that shouldn't be selected
#Returns:
#[ndarray int] med_ids: indices of the medoids selected
initMedoids <- function(X, n, init_type, exclude_ids){
  N<-nrow(X)
  D<-ncol(X)
  med_ids<-rep(-1,n)


  if(init_type=="kpp"){
    # Select a random point (not a boundary id)
    med_ids[1]<-sample(setdiff(1:nrows,boundary_ids),1)
    # Now fill the path with n randomly picked medoids
    boundary_coords<-input[boundary_ids,]
    for(i in 2:n){
      # Calculate all distances between points up to here and boundaries
      Xmed_dst<-pdist(input,rbind(input[med_ids[1:i],],input[boundary_ids,]))
      Xmed_dst<-as.matrix(Xmed_dst)^2
      D2<-apply(Xmed_dst,1,min)
      D2_n<-1.0/sum(D2)
    }


  } else if (init_type=='uniform'){
    stop("uniform type not implemented yet")
  } else {
    stop("init_type not recognized.")
  }

}

#Find the elbow in a fuction f, as the point on f with max distance from the line connecting f[0,:] and f[-1,:]
#[ndarray float] f: function(Nx2 array in the form [x, f(x)])
#[int] elb_id: index of the elbow
find_elbow <- function(f){

}
