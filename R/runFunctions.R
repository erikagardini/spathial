#' Select starting and ending points
#'
#' Get the coordinates of the starting and ending points
#'
#' @param X data points
#' @param X_labels labels of the data points
#' @param mode strategy for boundary selection
#' \itemize{
#'   \item 1 - centroids
#'   \item 2 - selected by the user
#'   \item 3 - insert the row name of the starting and ending points
#' }
#' @param from - starting class or row name of the starting point
#' @param to - ending class or row name of the ending point
#' @return list
#' \itemize{
#'   \item boundary ids - The indexes of the boundaries
#'   \item X - The new data matrix with the boundary
#'   \item X_labels - The new labels of the data matrix with the boundary labels
#' }
#' @export
spathial_boundary_ids <- function(X, X_labels, mode, from = NULL, to = NULL){
  if(mode == 2){
    X_2D <- spathial_2D(X)
    plot(X_2D$Y[,1],X_2D$Y[,2], pch=20,col="black",main="Click to select path start and end points")
    boundary_ids<-rownames(X)[identify(X,n=2,plot=FALSE)]
    points(
      X_2D$Y[boundary_ids,1], X_2D$Y[boundary_ids,2],pch="x",col="red",cex=4,
      xlab="Dimension 1",ylab="Dimension 2"
    )
  }else if(mode == 1){
    if(is.null(from) | is.null(to)){
      stop("You should insert the starting label and the ending label")
    }else{
      starting_centroid <- colMeans(X[which(X_labels == from),])
      ending_centroid <- colMeans(X[which(X_labels == to),])
      X <- rbind(X, starting_centroid, ending_centroid)
      rownames(X)[nrow(X):(nrow(X)-1)]<-c("Centroid2","Centroid1")
      X_labels <- c(X_labels, 0)
      X_labels <- c(X_labels, 0)
      names(X_labels)<-rownames(X)
      boundary_ids <- names(which(X_labels == 0,useNames = TRUE))
    }
  }else if(mode == 3){
    if(is.null(from) | is.null(to)){
      stop("You should insert the starting label and the ending label")
    }else if(!(from %in% rownames(X)) ){
      stop("from is not a valid class")
    }else if(!(to %in% rownames(X))){
      stop("to is not a valid class")
    }else{
      starting_point <- X[which(rownames(X) == from),]
      ending_point <- X[which(rownames(X) == to),]
      boundary_ids <- rownames(rbind(starting_point, ending_point))
    }
  }else{
    stop("Insert a valid mode")
  }

  outlist<-list(
    X=X,
    X_labels=X_labels,
    boundary_ids=boundary_ids
    )
  return(outlist)
}

#' Compute Principal Path
#'
#' Get the coordinates of the waypoints of the principal path
#'
#' @param X data points
#' @param boundary_ids starting and ending points
#' @param NC number of waypoints
#' @param prefiltering a boolean
#' @return ppath - spathial waypoints
#' @export
spathialWay <- function(X, boundary_ids, NC, prefiltering){
  if(prefiltering){
    ### Prefilter the data (function pp.rkm_prefilter)
    prefiltered<-rkm_prefilter(X,boundary_ids,plot_ax=TRUE)
    X<-prefiltered$X_filtered
    boundary_ids<-prefiltered$boundary_ids_filtered
    X_g<-prefiltered$X_garbage
    rm(prefiltered)
  }

  ### Initialize waypoints
  waypoint_ids<-initMedoids(X, NC, 'kpp', boundary_ids)
  waypoint_ids<-c(boundary_ids[1],waypoint_ids,boundary_ids[2])
  init_W<-X[waypoint_ids,]

  ### Annealing with rkm
  s_span<-pracma::logspace(5,-5)  #REMOVED ,n=NC -- the number of paths generated is different from the number of waypoint for each of them
  s_span<-c(s_span,0)
  #models<-array(data=NA,dim=c(length(s_span),NC+2,ncol(X)))
  #s<-s_span[1]

  models<-list()
  pb<-txtProgressBar(0,length(s_span),style=3)
  for(i in 1:length(s_span)){
    s<-s_span[i]
    W<-rkm(X,init_W,s,plot_ax=FALSE)
    init_W<-W
    models[[as.character(s)]]<-W
    #models[i,,]<-W
    setTxtProgressBar(pb,i)
  }
  W_dst_var <- rkm_MS_pathvar(models, s_span, X)
  s_elb_id <- find_elbow(cbind(s_span, W_dst_var))
  ppath <- models[[s_elb_id]]
  return(ppath)
}


#' Find labels
#'
#' Get the label of each waypoint accordin to the neighbourhood
#'
#' @param X data points
#' @param X_labels labels of the data points
#' @param ppath waypoints
#' @return ppath_labels - labels of the waypoints
#' @export
spathial_labels <- function(X, X_labels, ppath){
  library(class)
  X <- X[which(X_labels != 0),]
  X_labels <- X_labels[which(X_labels != 0)]
  ppath_no_centroids <- ppath[2:(nrow(ppath)-1), ]
  lbl <- knn(X, ppath_no_centroids, cl=X_labels, k=1)
  plot(c(1:length(lbl)), c(lbl), col=lbl, pch=19)
  return(lbl)
}

#' 2D spathial
#'
#' Get the 2D coordinates of each waypoint (using t-SNE algorithm for the dimensionality reduction)
#'
#' @param X data points
#' @param X_labels labels of the data points
#' @param ppath waypoints
#' @export
spathial_2D_plot <- function(X, X_labels, boundary_ids, ppath){

  ppath_labels <- array(data = -1, dim=(nrow(ppath)-2))

  total_labels <- c(X_labels, ppath_labels)
  all_points <- rbind(X, ppath[2:(nrow(ppath)-1), ])

  library(Rtsne())
  set.seed(1)
  tsne_res <- Rtsne(as.matrix(all_points), dims = 2, perplexity = 1)
  points_2D <- tsne_res$Y

  X_2D <- points_2D[which(total_labels != -1 & total_labels != 0),]
  boundary_ids_2D <- points_2D[which(rownames(X) == boundary_ids[1] | rownames(X) == boundary_ids[2]),]
  ppath_2D <- points_2D[which(total_labels == -1),]
  ppath_2D <- rbind(boundary_ids_2D[1,], ppath_2D, boundary_ids_2D[2,])

  plot(X_2D[,1],X_2D[,2], col=X_labels, pch=19)
  points(boundary_ids_2D[,1],boundary_ids_2D[,2], col="black", pch=3)
  lines(ppath_2D[,1], ppath_2D[,2],lwd=3,col="blue",type="o",pch=15)
}

#' Correlation
#'
#' ...
#'
#' @param ppath waypoints
#' @return corr - correlation along the path
#' @export
spathial_corr <- function(ppath){
  ppath_no_centroids <- ppath[2:(nrow(ppath)-1),]
  correlations <- apply(ppath, 2, function(x){
    cor(x, c(1:length(x)))
  })
  return(sort(correlations, decreasing = TRUE))
}
