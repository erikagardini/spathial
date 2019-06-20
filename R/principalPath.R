#' rkm
#'
#' Regularized K-means for principal path, MINIMIZER.
#'
#' @param X a data matrix
#' @param init_W initial waypoints matrix
#' @param s regularization parameter
#' @param plot_ax boolean, whether to plot a graph
#' @return W - final waypoints matrix
#' @references 'Finding Prinicpal Paths in Data Space', M.J.Ferrarotti, W.Rocchia, S.Decherchi, IEEE transactions on neural networks and learning systems 2018
#' @export
rkm <- function(X, init_W, s, plot_ax=FALSE){
  # Extract useful info from args
  N<-nrow(X)
  d<-ncol(X)
  NC<-nrow(init_W)-2

  # Construct boundary matrix
  boundary<-init_W[c(1,nrow(init_W)),]
  B<-matrix(0,nrow=NC,ncol=d)
  B[1,]<-as.numeric(boundary[1,])
  B[nrow(B),]<-as.numeric(boundary[2,])
  rownames(B)<-rownames(init_W)[2:(nrow(init_W)-1)]

  # Construct regularizer hessian
  diag1<-diag(NC)
  diag2<-diag(NC)
  diag(diag2)<-(-0.5)
  diag2<-cbind(rep(0,NC),diag2[,1:(ncol(diag2)-1)])
  diag3<-diag(NC)
  diag(diag3)<-(-0.5)
  diag3<-rbind(rep(0,NC),diag3[1:(nrow(diag3)-1),])
  AW<-diag1+diag2+diag3
  rownames(AW)<-colnames(AW)<-rownames(B)

  # Compute initial labels
  XW_dst<-pracma::distmat(
    as.matrix(X),
    as.matrix(init_W)
  )
  XW_dst<-XW_dst^2
  u<-colnames(XW_dst)[apply(XW_dst,1,which.min)]

  ### Iterate the minimizer
  converged<-FALSE
  it<-0
  while(!converged){
    it<-it+1
    #message("Iteration ",it)

    # Compute Cardinality
    W_card<-setNames(rep(0,NC+2),rownames(init_W))
    tu<-table(u)
    W_card[names(tu)]<-tu

    # Compute Centroid Matrix
    C<-matrix(NA,nrow=NC,ncol=d)
    rownames(C)<-rownames(B)
    colnames(C)<-colnames(X)
    for(i in 1:NC){ # NOTE: here I took some liberties, recheck (also it can be optimized and improved)
      ii<-rownames(C)[i]
      iiw<-which(u==ii)
      if(length(iiw)==1){
        centroid<-X[iiw,] # TODO I fixed this
      } else if(length(iiw)>1) {
        centroid<-apply(X[iiw,],2,sum)
      } else {
        centroid<-rep(0,ncol(X))
      }
      C[i,]<-centroid
    }

    # Construct K-means Heassian
    AX<-diag(W_card[2:(length(W_card)-1)])

    ### Update waypoints
    # Compute the (Moore-Penrose) pseudo-inverse of the multiplied Hessians
    pseudo<-MASS::ginv(AX+s*AW)
    # Multiply the centroid matrix with the s parameter with the boundaries, it's so clear it hurts
    csb<-C+0.5*s*B
    # Let's multiply these matrices, for the glory of the FaBiT
    W<-pseudo%*%csb
    rownames(W)<-rownames(B)
    colnames(W)<-colnames(X)
    W<-rbind(boundary[1,],W,boundary[2,])
    rownames(W)<-rownames(init_W)

    # Compute new labels
    XW_dst<-pracma::distmat(
      as.matrix(X),
      as.matrix(W)
    )
    XW_dst<-XW_dst^2
    u_new<-colnames(XW_dst)[apply(XW_dst,1,which.min)]

    # Check for convergence
    converged<-all(u_new==u)
    u<-u_new
    #converged<-(sum(u_new==u)<1) # NOTE: we can loosen the convergence requirements with >1
  }

  # Plot Progression
  #if(plot_ax){
   # plot(X[,57],X[,501],main=paste0("s=",s))
    #lines(W[,57], W[,501],lwd=3,col="red",type="o",pch=15)
  #}
  return(W)
}

#' rkm_prefilter
#'
#' Regularized K-means for principal path: prefiltering
#'
#' @param X a data matrix
#' @param boundary_ids names of the start and ending points, to be treated
#' separately
#' @param Nf parameter Nf (requires description)
#' @param k parameter k (requires description)
#' @param T parameter T (requires description)
#' @param plot_ax boolean, whether the post-filtering graph should be plotted
#' @return A list of objects
#' \itemize{
#'   \item X_filtered - The filtered data matrix
#'   \item boundary_ids_filtered - The filtered boundary ids
#'   \item X_garbage - The data matrix that was discarded
#' }
#' @export
rkm_prefilter <- function(X, boundary_ids, Nf=200, k=5, p=1000, T=0.1,
                          plot_ax=FALSE){
  N=nrow(X)
  d=ncol(X)
  # Pick Nf medoids with k-means++ and compute pairwise distance matrix
  med_ids<-initMedoids(X,n=Nf-2,init_type="kpp",boundary_ids=boundary_ids)
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
  medmed_dst_p[1,Nf]<-0 # NOTE: not sure why this is done
  medmed_dst_p[Nf,1]<-0 # NOTE: not sure why this is done

  # Find shortest path using dijkstra
  # NOTE: I didn't calculate the predecessors and simply got the shortest path
  g<-graph.adjacency(medmed_dst_p, weighted=TRUE)
  # NOTE: it is unnecessary to calculate the distances, but it could become
  # useful in the future, so I leave the code below:
  # path_dst<-igraph::distances(g,v=1,algorithm="dijkstra")
  path<-igraph::shortest_paths(g,from=1,to=ncol(medmed_dst_p))
  path<-names(path$vpath[[1]])

  # Filter out medoids too close to the shortest path
  T<-T*mean(medmed_dst)
  to_filter_ids<-c()
  for(i in path){
    to_filter_ids<-c(names(which(medmed_dst[i,]<T,useNames=TRUE)),to_filter_ids)
  }
  to_filter_ids<-unique(setdiff(to_filter_ids,path))
  to_keep_ids<-setdiff(rownames(medmed_dst),to_filter_ids)
  Xmed_dst<-pracma::distmat(
    as.matrix(X),
    as.matrix(X[to_keep_ids,])
  )
  Xmed_dst<-Xmed_dst^2
  rownames(Xmed_dst)<-rownames(X)
  colnames(Xmed_dst)<-to_keep_ids
  # Select minimum distance medoids
  u<-colnames(Xmed_dst)[apply(Xmed_dst,1,which.min)]
  filter_mask<-setNames(rep(FALSE,N),rownames(X))
  filter_mask[u%in%path]<-TRUE

  # NOTE: Convert boundary indices not necessary anymore,
  # I simply leveraged the object names power of R to avoid these tricks
  boundary_ids_filtered<-boundary_ids

  # Plot filter figure
  if(plot_ax){
    xcoord<-X[!filter_mask,1]
    ycoord<-X[!filter_mask,2]
    xrange<-abs(max(X[,1])-min(X[,1]))
    plot(xcoord,ycoord,
         xlab="Dimension 1",ylab="Dimension 2",
         main="Post-filtering Path plot",col="orange",pch=1,
         xlim=c(min(X[,1]),max(X[,1])+xrange/2)) # data filtered out
    points(X[filter_mask,1],X[filter_mask,2],pch=19,col="blue") # data kept
    points(X[med_ids,1],X[med_ids,2],col="red",pch=15) # filter medoids
    points(X[to_filter_ids, 1], X[to_filter_ids, 2],pch="x",cex=1) # filter medoids dropped
    lines(X[path,1], X[path,2],col="darkgreen",lwd=3,pch=19,type="o") # filter shortest path
    text(X[filter_mask,][boundary_ids_filtered,1],
         X[filter_mask,][boundary_ids_filtered,2],
         labels=c("B1","B2"),font=2,cex=1.3,col="grey") # Boundary samples
    legend("right",
           legend=c(
             "data filtered out",
             "data kept",
             "filter medoids",
             "filter medoids dropped",
             "filter shortest path",
             "boundary samples"
           ),bg="#FFFFFF22",
           pch=c(1,19,15,4,19,66),
           col=c("orange","blue","red","black","darkgreen","grey"),
           lty=c(0,0,0,0,1,0),lwd=c(0,0,0,0,2,0)
    )
  }
  # Return ouput
  outlist<-list(
    X_filtered=X[filter_mask,],
    boundary_ids_filtered=boundary_ids_filtered,
    X_garbage=X[!filter_mask,]
  )
  return(outlist)
}

#' Model Selection
#'
#' Regularized K-menans for principal path, MODEL SELECTION, variance on waypoints interdistance.
#'
#' @param modelsmatrix with path models, shape, N_models x N x (NC+2)
#' @param s_span array with values of the reg parameter for each model (sorted in decreasing order, with 0 as last value)
#' @param X data matrix
#' @return W_dst_var: array with values of variance for each model
#' @export
rkm_MS_pathvar <- function(models, s_span, X){
  W_dst_var <- array(data = 0.0, dim = length(models))
  for(i in (1:length(models))){
    W = models[[i]]
    W=as.matrix(W)
    res <- W[2:dim(W)[1],] - W[1:dim(W)[1]-1,]
    W_dst=apply(res, 1, function(x) norm(as.matrix(x),"F"))
    W_dst_var[i] = var(W_dst)
  }
  return(W_dst_var)
}
