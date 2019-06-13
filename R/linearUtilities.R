#Initialize NC medoids with init_type rational.
#[ndarray float] X: data matrix
#[int] n: number of medoids to be selected
#[string] init_type: rational to be used ('uniform' = randomly selected with uniform distribution, 'kpp' = k-means++ algorithm)
#[ndarray int] exclude_ids: blacklisted ids that shouldn't be selected
#[ndarray int] med_ids: indices of the medoids selected
initMedoids <- function(X, n, init_type, exclude_ids){

}

#Find the elbow in a fuction f, as the point on f with max distance from the line connecting f[0,:] and f[-1,:]
#[ndarray float] f: function(Nx2 array in the form [x, f(x)])
#[int] elb_id: index of the elbow
find_elbow <- function(f){

}
