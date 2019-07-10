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
#' @param from starting class or row name of the starting point
#' @param to ending class or row name of the ending point
#' @return A list of objects
#' \itemize{
#'   \item boundary ids: the indexes of the boundaries
#'   \item X: the new data matrix with the boundary
#'   \item X_labels: the new labels of the data matrix with the boundary labels
#' }
#' @export
spathialBoundaryIds <- function(X, X_labels, mode, from = NULL, to = NULL){
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
    }else if(!(from %in% rownames(X)) ){
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
#' @param X a data matrix
#' @param X_labels labels of the data points
#' @param boundary_ids names of the start and ending points, to be treated separately
#' @return A list of objects
#' \itemize{
#'   \item X_filtered: the filtered data matrix
#'   \item X_labels_filtered: the filtered labels of the data points
#'   \item boundary_ids_filtered: the filtered boundary ids
#'   \item X_garbage: the data matrix that was discarded
#'   \item X_garbage_filtered: the filterad labels that was discarded
#' }
#' @export
spathialPrefiltering <- function(X, X_labels, boundary_ids){
  prefiltered<-rkm_prefilter(X, boundary_ids)
  X_filtered<-X[prefiltered$filter_mask,]
  X_labels_filtered <- X_labels[prefiltered$filter_mask]
  X_garbage<-X[!prefiltered$filter_mask,]
  X_labels_garbage<-X_labels[!prefiltered$filter_mask]
  X_labels_garbage <- sapply(X_labels_garbage, function(x){
    return(-2)
  })
  boundary_ids<-prefiltered$boundary_ids_filtered
  rm(prefiltered)

  outlist<-list(
    X_filtered=X_filtered,
    X_labels_filtered=X_labels_filtered,
    X_garbage=X_garbage,
    X_labels_garbage=X_labels_garbage,
    boundary_ids=boundary_ids
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
#' @param neighbors the number of desired nearest neighbours
#' @return A list of objects
#' \itemize{
#'   \item ppath: spathial waypoints
#'   \item ppath_parturbed: all the spathial waypoints for each perturbation
#'}
#' @export
spathialWay <- function(X, X_labels, boundary_ids, NC, neighbors = NULL){
  if(is.null(neighbors)){
    neighbors <- 1
  }
  if(neighbors == 1){
    ppath <- compute_spathial(X, boundary_ids, NC)
    colnames(ppath) <- colnames(X)
    perturbed_path <- NULL
  }
  else{
    starting_class <- X_labels[which(rownames(X) == boundary_ids[1])]
    ending_class <- X_labels[which(rownames(X) == boundary_ids[2])]

    element_starting_class <- X[which(X_labels == starting_class),]
    element_ending_class <- X[which(X_labels == ending_class),]

    starting_class_neighbour <- find_nearest_points(X[which(rownames(X) == boundary_ids[1]),], element_starting_class, neighbors)
    ending_class_neighbour <- find_nearest_points(X[which(rownames(X) == boundary_ids[2]),], element_ending_class, neighbors)

    perturbed_path <- lapply(starting_class_neighbour, function(x){
      lapply(ending_class_neighbour, function(y){
        boundary_ids <- c(x, y)
        perturbed <- compute_spathial(X, boundary_ids, NC)
        colnames(perturbed) <- colnames(X)
        return(perturbed)
      })
    })
    ppath <- perturbed_path[[1]][[1]]
    colnames(ppath) <- colnames(X)
  }

  outlist<-list(
    ppath=ppath,
    perturbed_path=perturbed_path
  )
  return(outlist)
}

#' Find labels
#'
#' Get the label of each waypoint accordin to the neighbourhood
#'
#' @param X data points
#' @param X_labels labels of the data points
#' @param spathial_res A list of objects
#' \itemize{
#'   \item ppath: principal path from the starting point to the ending point
#'   \item ppath_perturbed: all the perturbed paths
#' }
#' @return ppath_labels - labels of the waypoints
#' @export
spathialLabels <- function(X, X_labels, spathial_res){
  ppath <- spathial_res$ppath
  X_labels <- X_labels[which(! grepl("Centroid", rownames(X)))]
  X <- X[which(! grepl("Centroid", rownames(X))),]
  ppath_no_centroids <- ppath[2:(nrow(ppath)-1), ]
  lbl <- class::knn(X, ppath_no_centroids, cl=X_labels, k=1)
  plot(c(1:length(lbl)), c(lbl), col=lbl, pch=19)
  return(lbl)
}

#' 2D spathial
#'
#' Get the 2D coordinates of each waypoint (using t-SNE algorithm for the dimensionality reduction)
#'
#' @param X data points
#' @param X_labels labels of the data points
#' @param boundary_ids waypoints
#' @param spathial_res A list of objects
#' \itemize{
#'   \item ppath: principal path from the starting point to the ending point
#'   \item ppath_perturbed: all the perturbed paths
#' }
#' @param perplexity_value the value for TSNE perplexity (default is nrsamples*3/50)
#' @param X_garbage the data points removed during the filtering (if the prefiltering is done)
#' @param X_labels_garbage the labels of the data points removed during the filtering (if the prefiltering is done)
#' @export
spathialPlot2D <- function(X, X_labels, boundary_ids, spathial_res, perplexity_value=NULL, X_garbage=NULL, X_labels_garbage=NULL){
  set.seed(1)

  # Exception handler if the user doesn't specificy perplexity
  if(is.null(perplexity_value)){
    perplexity_value<-ceiling(nrow(X)*3/50)
    #message("Perplexity is ",perplexity_value)
  }
  ppath <- spathial_res$ppath
  ppath <- ppath[2:(nrow(ppath)-1),]
  rownames(ppath) <- paste("ppath",1:nrow(ppath))
  ppath_labels <- array(data = -1, dim=(nrow(ppath)))
  total_labels <- c(X_labels, ppath_labels)
  all_points <- rbind(X, ppath)
  if(!is.null(X_garbage)){
    total_labels <- c(total_labels, X_labels_garbage)
    all_points <- rbind(all_points, X_garbage)
  }

  tsne_res <- Rtsne::Rtsne(as.matrix(all_points), dims = 2, perplexity = perplexity_value)
  points_2D <- tsne_res$Y

  X_2D <- points_2D[which(total_labels != -1 & total_labels != 0 & total_labels != -2),]
  boundary_ids_2D <- points_2D[which(rownames(X) == boundary_ids[1] | rownames(X) == boundary_ids[2]),]
  ppath_2D <- points_2D[which(total_labels == -1),]
  ppath_2D <- rbind(boundary_ids_2D[1,], ppath_2D, boundary_ids_2D[2,])
  if(!is.null(X_garbage)){
    X_garbage_2D <- points_2D[which(total_labels == -2),]
  }

  plot(X_2D[,1],X_2D[,2], col=X_labels, pch=19)
  points(boundary_ids_2D[,1],boundary_ids_2D[,2], col="black", pch=3)
  lines(ppath_2D[,1], ppath_2D[,2],lwd=3,col="blue",type="o",pch=15)
  if(!is.null(X_garbage)){
    points(X_garbage_2D[,1],X_garbage_2D[,2], col="gray", pch=3)
  }
}

#' Correlation
#'
#' Get how much the features correlate with the path
#'
#' @param spathial_res A list of objects
#' \itemize{
#'   \item ppath: principal path from the starting point to the ending point
#'   \item ppath_perturbed: all the perturbed paths
#' }
#' @return A list of objects
#' \itemize{
#'   \item ppath: spathial waypoints
#'   \item ppath_parturbed: all the spathial waypoints for each perturbation
#'}
#' @export
spathialCorrelation <- function(spathial_res){
  #Mean of z scores (with Fisher transformation)
  if(!is.null(spathial_res$perturbed_path)){
    z_scores_perturbed_path <- lapply(spathial_res$perturbed_path, function(x){
      lapply(x, function(y){
        z_scores <- apply(y, 2, function(z){
          DescTools::FisherZ(cor(z, c(1:length(z))))
        })
      })
    })

    count <- 0
    sum <- array(data=0, dim=c(ncol(spathial_res$ppath)))
    for(i in (1:length(z_scores_perturbed_path))){
      A <- z_scores_perturbed_path[[i]]
      for(j in (1:length(A))){
        B <- A[[j]]
        for(k in (1:length(B))){
          sum[k] <- sum[k] +  B[k]
        }
        count <- count + 1
      }
    }
    fisher <- sapply(sum, function(x){
      z_avg <- x/count
      return(DescTools::FisherZInv(z_avg))
    })
    fisher <- as.list(fisher)
    names(fisher) <- colnames(spathial_res$ppath)

    #Simple mean of correlation
    # correlations <- lapply(spathial_res$perturbed_path, function(x){
    #   lapply(x, function(y){
    #     corr <- apply(y, 2, function(z){
    #       cor(z, c(1:length(z)))
    #     })
    #   })
    # })
    # correlations <- colMeans(correlations)
    correlations <- NULL

  }else{
    correlations <- lapply(spathial_res$ppath, function(x){
      cor(x, c(1:length(x)))
    })
    names(correlations) <- colnames(spathial_res$ppath)
    fisher <- NULL
  }

  outlist<-list(
    correlations=correlations,
    fisher=fisher
  )
  return(outlist)
}
