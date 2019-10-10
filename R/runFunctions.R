#' Select starting and ending points
#'
#' Get the coordinates of the starting and ending points
#'
#' @param X data points
#' @param X_labels labels of the data points
#' @param mode strategy for boundary selection
#' \itemize{
#'   \item 1 - selected by the user
#'   \item 2 - centroids
#'   \item 3 - insert the row name of the starting and ending points
#' }
#' @param from starting class or row name of the starting point
#' @param to ending class or row name of the ending point
#' @return A list of objects
#' \itemize{
#'   \item boundary ids: the indexes of the boundaries
#'   \item X: the new data matrix with the boundary
#'   \item X_labels: the new labels of the data matrix with the boundary labels
#' }
#' @examples
#' # Load data matrix X
#' load(system.file('extdata','X.rda',package='spathial',mustWork=TRUE))
#' # Load description vector X_labels
#' load(system.file('extdata','X_labels.rda',package='spathial',mustWork=TRUE))
#' # Run spathialBoundary
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' @export
spathialBoundaryIds <- function(X, X_labels, mode = 1, from = NULL, to = NULL){
  if(mode == 1){
    if(ncol(X) < 2){
      stop("X should have at least 2 columns")
    }else if(ncol(X) == 2){
      colors <- rainbow(length(table(X_labels)))
      colors_labels <- sapply(X_labels, function(x){colors[x]})

      plot(X[,1],X[,2], col=colors_labels, pch=as.character(X_labels), main="Click to select path start and end points")
      boundary_ids<-rownames(X)[identify(X,n=2,plot=FALSE)]
      points(
        X[which(rownames(X) == boundary_ids[1]),1], X[which(rownames(X) == boundary_ids[1]),2],pch="x",col="black",cex=4,
        xlab="Dimension 1",ylab="Dimension 2"
      )
      points(
        X[which(rownames(X) == boundary_ids[2]),1], X[which(rownames(X) == boundary_ids[2]),2],pch="x",col="black",cex=4,
        xlab="Dimension 1",ylab="Dimension 2"
      )
    }else{
      tsne_res <- Rtsne::Rtsne(X, dims = 2, perplexity = 30)
      X_2D <- tsne_res$Y

      colors <- rainbow(length(table(X_labels)))
      colors_labels <- sapply(X_labels, function(x){colors[x]})

      plot(X_2D[,1],X_2D[,2], col=colors_labels, pch=as.character(X_labels), main="Click to select path start and end points")
      boundary_ids<-rownames(X)[identify(X_2D,n=2,plot=FALSE)]
      points(
        X_2D[which(rownames(X) == boundary_ids[1]),1], X_2D[which(rownames(X) == boundary_ids[1]),2],pch="x",col="black",cex=4,
        xlab="Dimension 1",ylab="Dimension 2"
      )
      points(
        X_2D[which(rownames(X) == boundary_ids[2]),1], X_2D[which(rownames(X) == boundary_ids[2]),2],pch="x",col="black",cex=4,
        xlab="Dimension 1",ylab="Dimension 2"
      )
    }
  }else if(mode == 2){
    if(is.null(from) | is.null(to)){
      stop("You should insert the starting label and the ending label")
    }else if(!(from %in% X_labels)){
      stop("from is not a valid class")
    }else if(!(to %in% X_labels)){
      stop("to is not a valid class")
    }else{
      starting_centroid <- colMeans(X[which(X_labels == from),], na.rm = TRUE)
      ending_centroid <- colMeans(X[which(X_labels == to),], na.rm = TRUE)
      X <- rbind(X, starting_centroid, ending_centroid)
      rownames(X)[nrow(X):(nrow(X)-1)]<-c("Centroid2","Centroid1")
      X_labels <- c(X_labels, from)
      X_labels <- c(X_labels, to)
      names(X_labels)<-rownames(X)
      boundary_ids <- rownames(X[grep("Centroid", rownames(X)),])
    }
  }else if(mode == 3){
    if(is.null(from) | is.null(to)){
      stop("You should insert the starting label and the ending label")
    }else if(!(from %in% rownames(X))){
      stop("from is not an existing sample")
    }else if(!(to %in% rownames(X))){
      stop("to is not an existing sample")
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

#' Prefilter data
#'
#' Regularized K-means for principal path: prefiltering
#'
#' @param X data points
#' @param X_labels labels of the data points
#' @param boundary_ids names of the start and ending points, to be treated separately
#' @return A list of objects
#' \itemize{
#'   \item mask: indexes of the data points to preserv
#'   \item boundary_ids: the filtered boundary ids
#' }
#' @examples
#' # Load data matrix X
#' load(system.file('extdata','X.rda',package='spathial',mustWork=TRUE))
#' # Load description vector X_labels
#' load(system.file('extdata','X_labels.rda',package='spathial',mustWork=TRUE))
#' # Run spathialBoundary
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' # Run spathial spathialPrefilterinh with the output of the function spathialBoundaryIds
#' filterRes <- spathialPrefiltering(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids)
#' @export
spathialPrefiltering <- function(X, X_labels, boundary_ids){
  prefiltered<-rkm_prefilter(X, boundary_ids)

  outlist<-list(
    mask=prefiltered$filter_mask,
    boundary_ids=prefiltered$boundary_ids_filtered
  )
  return(outlist)
}

#' Compute Principal Path
#'
#' Get the coordinates of the waypoints of the principal path
#'
#' @param X data points
#' @param X_labels labels of the data points
#' @param boundary_ids starting and ending points
#' @param NC number of waypoints
#' @return spathial_res: spathial waypoints
#' @examples
#' #EXAMPLE 1
#' # Load data matrix X
#' load(system.file('extdata','X.rda',package='spathial',mustWork=TRUE))
#' # Load description vector X_labels
#' load(system.file('extdata','X_labels.rda',package='spathial',mustWork=TRUE))
#' # Run spathialBoundary
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' #Set the number of waypoints
#' NC <- 20
#' # Run spathialWay without prefilterint
#' spathial_res <- spathialWay(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids, NC)
#' spathial_res
#'
#' #EXAMPLE 2
#' # Load data matrix X
#' load(system.file('extdata','X.rda',package='spathial',mustWork=TRUE))
#' # Load description vector X_labels
#' load(system.file('extdata','X_labels.rda',package='spathial',mustWork=TRUE))
#' # Run spathialBoundary
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' # Run spathial spathialPrefilterinh with the output of the function spathialBoundaryIds
#' filterRes <- spathialPrefiltering(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids)
#' #Set the number of waypoints
#' NC <- 20
#' # Run spathialWay with the filtered data
#' spathial_res <- spathialWay(boundaryRes$X[filterRes$mask,], boundaryRes$X_labels[filterRes$mask], filterRes$boundary_ids, NC)
#' spathial_res
#' @export
spathialWay <- function(X, X_labels, boundary_ids, NC=50){
  spathial_res <- compute_spathial(X, boundary_ids, NC)
  colnames(spathial_res) <- colnames(X)

  return(spathial_res)
}

#' Find labels
#'
#' Get the label of each waypoint accordin to the neighbourhood
#'
#' @param X data points
#' @param X_labels labels of the data points
#' @param spathial_res: principal path from the starting point to the ending point
#' @return ppath_labels: labels of the waypoints
#' @examples
#' #EXAMPLE 1
#' # Load data matrix X
#' load(system.file('extdata','X.rda',package='spathial',mustWork=TRUE))
#' # Load description vector X_labels
#' load(system.file('extdata','X_labels.rda',package='spathial',mustWork=TRUE))
#' # Run spathialBoundary
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' #Set the number of waypoints
#' NC <- 20
#' # Run spathialWay without prefiltering
#' spathial_res <- spathialWay(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids, NC)
#' #Run spathialLabels with spathial_res
#' labels <- spathialLabels(boundaryRes$X, boundaryRes$X_labels, spathial_res)
#' labels
#'
#' #EXAMPLE 2
#' # Load data matrix X
#' load(system.file('extdata','X.rda',package='spathial',mustWork=TRUE))
#' # Load description vector X_labels
#' load(system.file('extdata','X_labels.rda',package='spathial',mustWork=TRUE))
#' # Run spathialBoundary
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' # Run spathial spathialPrefilterinh with the output of the function spathialBoundaryIds
#' filterRes <- spathialPrefiltering(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids)
#' #Set the number of waypoints
#' NC <- 20
#' # Run spathialWay with the filtered data
#' spathial_res <- spathialWay(boundaryRes$X[filterRes$mask,], boundaryRes$X_labels[filterRes$mask], filterRes$boundary_ids, NC)
#' #Run spathialLabels with spathial_res
#' labels <- spathialLabels(boundaryRes$X[filterRes$mask,], boundaryRes$X_labels[filterRes$mask], spathial_res)
#' labels
#' @export
spathialLabels <- function(X, X_labels, spathial_res){
  X_labels <- X_labels[which(! grepl("Centroid", rownames(X)))]
  X <- X[which(! grepl("Centroid", rownames(X))),]
  ppath_no_centroids <- spathial_res[2:(nrow(spathial_res)-1), ]
  ppath_labels <- class::knn(X, ppath_no_centroids, cl=X_labels, k=1)
  return(ppath_labels)
}

#' 2D spathial
#'
#' Get the 2D coordinates of each waypoint (using t-SNE algorithm for the dimensionality reduction)
#'
#' @param X data points
#' @param X_labels labels of the data points
#' @param boundary_ids waypoints
#' @param spathial_res: principal path from the starting point to the ending point
#' @param perplexity_value the value for TSNE perplexity (default is nrsamples*3/50)
#' @param mask the mask of the sample to preserve (when prefiltering is computed)
#' @param ... Parameters which will be inherited by plot()
#' @examples
# Load data matrix X
#' load(system.file('extdata','X.rda',package='spathial',mustWork=TRUE))
#' # Load description vector X_labels
#' load(system.file('extdata','X_labels.rda',package='spathial',mustWork=TRUE))
#' # Run spathialBoundary
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' #Set the number of waypoints
#' NC <- 20
#' # Run spathialWay without prefiltering
#' spathial_res <- spathialWay(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids, NC)
#' #Run spathialPlot with spathial_res
#' spathialPlot(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids, spathial_res, perplexity_value=30)
#'
#' #EXAMPLE 2
#' # Load data matrix X
#' load(system.file('extdata','X.rda',package='spathial',mustWork=TRUE))
#' # Load description vector X_labels
#' load(system.file('extdata','X_labels.rda',package='spathial',mustWork=TRUE))
#' # Run spathialBoundary
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' # Run spathial spathialPrefilterinh with the output of the function spathialBoundaryIds
#' filterRes <- spathialPrefiltering(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids)
#' #Set the number of waypoints
#' NC <- 20
#' # Run spathialWay with the filtered data
#' spathial_res <- spathialWay(boundaryRes$X[filterRes$mask,], boundaryRes$X_labels[filterRes$mask], filterRes$boundary_ids, NC)
#' #Run spathialPlot with spathial_res
#' spathialPlot(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids, spathial_res, perplexity_value=30, filterRes$mask)
#' @export
spathialPlot <- function(X, X_labels, boundary_ids, spathial_res, perplexity_value=NULL, mask=NULL, ...){
  set.seed(1)
  if(ncol(X) == 2){
    colors <- rainbow(length(table(X_labels)))
    colors_labels <- sapply(X_labels, function(x){colors[x]})
    boundaries <- X[which(rownames(X) == boundary_ids[1] | rownames(X) == boundary_ids[2]),]

    plot(X[,1],X[,2], col=colors_labels, xlab="tsne1", ylab="tsne2", ypch=as.character(X_labels),...)
    points(boundaries[,1],boundaries[,2], pch="x",col="black",cex=4)
    lines(spathial_res[,1], spathial_res[,2],lwd=3,col="red",type="o",pch=15)
    if(!is.null(mask)){
      X_garbage <- X[!mask,]
      X_labels_garbage <- X_labels[!mask]
      points(X_garbage[,1],X_garbage[,2], col="gray", pch="x")
    }
  }else{
    if(is.null(perplexity_value)){
      perplexity_value<-ceiling(nrow(X)*3/50)
      #message("Perplexity is ",perplexity_value)
    }
    rownames(spathial_res) <- paste("ppath",1:nrow(spathial_res))
    ppath_labels <- array(data = -1, dim=(nrow(spathial_res)))
    total_labels <- c(X_labels, ppath_labels)
    all_points <- rbind(X, spathial_res)

    tsne_res <- Rtsne::Rtsne(as.matrix(all_points), dims = 2, perplexity = perplexity_value, check_duplicates=FALSE)
    points_2D <- tsne_res$Y

    boundary_ids_2D <- points_2D[which(rownames(X) == boundary_ids[1] | rownames(X) == boundary_ids[2]),]
    ppath_2D <- points_2D[which(total_labels == -1),]
    #ppath_2D <- rbind(boundary_ids_2D[1,], ppath_2D, boundary_ids_2D[2,])

    points_2D <- points_2D[which(total_labels != -1),]

    if(!is.null(mask)){
      X_2D <- points_2D[mask,]
      X_garbage_2D <- points_2D[!mask, ]
      X_labels <- X_labels[mask]
    }

    colors <- rainbow(length(table(X_labels)))
    colors_labels <- sapply(X_labels, function(x){colors[x]})

    plot(points_2D[,1],points_2D[,2], xlab="tsne1", ylab="tsne2", col=colors_labels, pch=as.character(X_labels),...)
    points(boundary_ids_2D[,1],boundary_ids_2D[,2], pch="x",col="black",cex=4)
    lines(ppath_2D[,1], ppath_2D[,2],lwd=3,col="blue",type="o",pch=15)
    if(!is.null(mask)){
      points(X_garbage_2D[,1],X_garbage_2D[,2], col="gray", pch="x")
    }
  }
}

#' Correlation
#'
#' Get how much the features correlate with the path
#'
#' @param spathial_res: principal path from the starting point to the ending point
#' @return A list of objects
#' \itemize{
#'   \item correlations: Pearson's correlation coefficients between ea
#'   ch feature and the path (when ppath_perturbed is not NULL, a Fisher-integrated correlation coefficient is provided)
#'   \item ranks: ranks of associations between the n features and the path (when ppath_perturbed is not NULL, the mean of the ranks is provided)
#'}
#' @examples
#' # Load data matrix X
#' load(system.file('extdata','X.rda',package='spathial',mustWork=TRUE))
#' # Load description vector X_labels
#' load(system.file('extdata','X_labels.rda',package='spathial',mustWork=TRUE))
#' # Run spathialBoundary
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' # Run spathial spathialPrefilterinh with the output of the function spathialBoundaryIds
#' filterRes <- spathialPrefiltering(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids)
#' #Set the number of waypoints
#' NC <- 20
#' # Run spathialWay with the filtered data
#' spathial_res <- spathialWay(boundaryRes$X[filterRes$mask,], boundaryRes$X_labels[filterRes$mask], filterRes$boundary_ids, NC)
#' #Run spathialStatistics with spathial_res
#' statistics <- spathialStatistics(spathial_res)
#' @export
spathialStatistics <- function(spathial_res){
  correlations <- sapply(spathial_res, function(x){
    if(sd(x) == 0){
      return(0)
    }else{
      cor(x, c(1:length(x)))
    }
  })
  correlations<-unlist(correlations)
  names(correlations)<-colnames(spathial_res)
  ranks <- rank(-correlations)
  names(ranks)<-colnames(spathial_res)

  outlist<-list(
    correlations=correlations,
    ranks=ranks
  )
  return(outlist)
}

