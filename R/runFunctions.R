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
#' load("~/spathial/inst/extdata/X.rda")
#' load("~/spathial/inst/extdata/X_labels.rda")
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' @export
spathialBoundaryIds <- function(X, X_labels, mode = 1, from = NULL, to = NULL){
  if(mode == 1){
    if(ncol(X) < 2){
      stop("X should have at least 2 columns")
    }else if(ncol(X) == 2){
      colors <- rainbow(length(table(X_labels)))
      colors_labels <- sapply(X_labels, function(x){colors[x]})

      plot(X[,1],X[,2], col=color_labels, pch=as.character(X_labels), main="Click to select path start and end points")
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
#' @examples
#' load("~/spathial/inst/extdata/X.rda")
#' load("~/spathial/inst/extdata/X_labels.rda")
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
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
#' @param neighbors the number of desired nearest neighbours
#' @return A list of objects
#' \itemize{
#'   \item ppath: spathial waypoints
#'   \item ppath_parturbed: all the spathial waypoints for each perturbation
#'}
#' @examples
#' load("~/spathial/inst/extdata/X.rda")
#' load("~/spathial/inst/extdata/X_labels.rda")
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' NC <- 20
#' spathialRes <- spathialWay(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids, NC, neighbors = 1)
#'
#' load("~/spathial/inst/extdata/X.rda")
#' load("~/spathial/inst/extdata/X_labels.rda")
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' filterRes <- spathialPrefiltering(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids)
#' NC <- 20
#' spathialRes <- spathialWay(boundaryRes$X[filterRes$mask,], boundaryRes$X_labels[filterRes$mask], filterRes$boundary_ids, NC, neighbors = 1)
#' @export
spathialWay <- function(X, X_labels, boundary_ids, NC=50, neighbors = 0){
  if(neighbors == 0){
    ppath <- compute_spathial(X, boundary_ids, NC)
    colnames(ppath) <- colnames(X)
    perturbed_paths <- NULL
  }
  else{
    starting_class <- X_labels[which(rownames(X) == boundary_ids[1])]
    ending_class <- X_labels[which(rownames(X) == boundary_ids[2])]

    element_starting_class <- X[which(X_labels == starting_class),]
    element_ending_class <- X[which(X_labels == ending_class),]

    starting_class_neighbour <- find_nearest_points(X[which(rownames(X) == boundary_ids[1]),], element_starting_class, neighbors+1)
    ending_class_neighbour <- find_nearest_points(X[which(rownames(X) == boundary_ids[2]),], element_ending_class, neighbors+1)

    tot_path <- length(starting_class_neighbour) * length(ending_class_neighbour)
    i <- 1
    perturbed_paths <- lapply(starting_class_neighbour, function(x){
      lapply(ending_class_neighbour, function(y){
        message(paste0("Path ", i, "/", tot_path))
        i <<- i + 1
        boundary_ids <- c(x, y)
        perturbed <- compute_spathial(X, boundary_ids, NC)
        colnames(perturbed) <- colnames(X)
        return(perturbed)
      })
    })
    ppath <- perturbed_paths[[1]][[1]]
    colnames(ppath) <- colnames(X)
  }

  outlist<-list(
    ppath=ppath,
    perturbed_paths=perturbed_paths
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
#' @examples
#' load("~/spathial/inst/extdata/X.rda")
#' load("~/spathial/inst/extdata/X_labels.rda")
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' NC <- 20
#' spathialRes <- spathialWay(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids, NC, neighbors = 1)
#' labels <- spathialLabels(boundaryRes$X, boundaryRes$X_labels, spathialRes)
#'
#' load("~/spathial/inst/extdata/X.rda")
#' load("~/spathial/inst/extdata/X_labels.rda")
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' filterRes <- spathialPrefiltering(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids)
#' NC <- 20
#' spathialRes <- spathialWay(boundaryRes$X[filterRes$mask,], boundaryRes$X_labels[filterRes$mask], filterRes$boundary_ids, NC, neighbors = 1)
#' labels <- spathialLabels(boundaryRes$X[filterRes$mask,], boundaryRes$X_labels[filterRes$mask], spathialRes)
#' @export
spathialLabels <- function(X, X_labels, spathial_res){
  ppath <- spathial_res$ppath
  X_labels <- X_labels[which(! grepl("Centroid", rownames(X)))]
  X <- X[which(! grepl("Centroid", rownames(X))),]
  ppath_no_centroids <- ppath[2:(nrow(ppath)-1), ]
  lbl <- class::knn(X, ppath_no_centroids, cl=X_labels, k=1)
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
#' @param ... Parameters which will be inherited by plot()
#' \itemize{
#'   \item ppath: principal path from the starting point to the ending point
#'   \item ppath_perturbed: all the perturbed paths
#' }
#' @param perplexity_value the value for TSNE perplexity (default is nrsamples*3/50)
#' @param mask the mask of the sample to preserve (when prefiltering is computed)
#' @examples
#' load("~/spathial/inst/extdata/X.rda")
#' load("~/spathial/inst/extdata/X_labels.rda")
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' NC <- 20
#' spathialRes <- spathialWay(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids, NC, neighbors = 1)
#' spathialPlot(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids, spathialRes, perplexity_value=30)
#'
#' load("~/spathial/inst/extdata/X.rda")
#' load("~/spathial/inst/extdata/X_labels.rda")
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' filterRes <- spathialPrefiltering(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids)
#' NC <- 20
#' spathialRes <- spathialWay(boundaryRes$X[filterRes$mask,], boundaryRes$X_labels[filterRes$mask], filterRes$boundary_ids, NC, neighbors = 1)
#' spathialPlot(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids, spathialRes, perplexity_value=30, filterRes$mask)
#' @export
spathialPlot <- function(X, X_labels, boundary_ids, spathial_res,
                         perplexity_value=NULL, mask=NULL, ...){
  set.seed(1)
  if(is.null(spathial_res$perturbed_paths)){
    if(ncol(X) == 2){
      ppath <- spathial_res$ppath

      colors <- rainbow(length(table(X_labels)))
      colors_labels <- sapply(X_labels, function(x){colors[x]})
      boundaries <- X[which(rownames(X) == boundary_ids[1] | rownames(X) == boundary_ids[2]),]

      plot(X[,1],X[,2], col=colors_labels, pch=as.character(X_labels),...)
      points(boundaries[,1],boundaries[,2], pch="x",col="black",cex=4)
      lines(ppath[,1], ppath[,2],lwd=3,col="red",type="o",pch=15)
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
      ppath <- spathial_res$ppath
      #ppath <- ppath[2:(nrow(ppath)-1),]
      rownames(ppath) <- paste("ppath",1:nrow(ppath))
      ppath_labels <- array(data = -1, dim=(nrow(ppath)))
      total_labels <- c(X_labels, ppath_labels)
      all_points <- rbind(X, ppath)

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

      plot(points_2D[,1],points_2D[,2], col=colors_labels, pch=as.character(X_labels),...)
      points(boundary_ids_2D[,1],boundary_ids_2D[,2], pch="x",col="black",cex=4)
      lines(ppath_2D[,1], ppath_2D[,2],lwd=3,col="blue",type="o",pch=15)
      if(!is.null(mask)){
        points(X_garbage_2D[,1],X_garbage_2D[,2], col="gray", pch="x")
      }
    }
  }
  else{
    if(ncol(X) == 2){
      ppaths <- spathial_res$perturbed_paths

      colors <- rainbow(length(table(X_labels)))
      colors_labels <- sapply(X_labels, function(x){colors[x]})
      boundaries <- X[which(rownames(X) == boundary_ids[1] | rownames(X) == boundary_ids[2]),]

      plot(X[,1],X[,2], col=colors_labels, pch=as.character(X_labels),...)
      points(boundaries[,1],boundaries[,2], pch="x",col="black",cex=4)

      colors <- rainbow(length(ppaths) + length(ppaths[[1]]))
      i = 1
      for(el in ppaths){
        for(path in el){
          lines(path[,1], path[,2],lwd=3,col=colors[i],type="o",pch=i)
          i <- i + 1
        }
      }

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
      perturbed_paths <- spathial_res$perturbed_paths

      perturbed_paths_matrix <- matrix(nrow = 0, ncol = ncol(X))
      perturbed_path_labels <- matrix(nrow=0, ncol=1)

      i <- -1
      for(el in perturbed_paths){
        for(path in el){
          path <- as.matrix(path)
          perturbed_paths_matrix <- rbind(perturbed_paths_matrix, path)
          perturbed_path_labels <- c(perturbed_path_labels, rep(i, nrow(path)))
          i <- i - 1
        }
      }

      total_labels <- c(X_labels, perturbed_path_labels)
      all_points <- rbind(X, perturbed_paths_matrix)

      tsne_res <- Rtsne::Rtsne(as.matrix(all_points), dims = 2, perplexity = perplexity_value, check_duplicates = FALSE)
      points_2D <- tsne_res$Y

      boundary_ids_2D <- points_2D[which(rownames(X) == boundary_ids[1] | rownames(X) == boundary_ids[2]),]
      ppath_2D <- points_2D[which(total_labels < 0),]
      points_2D <- points_2D[which(total_labels > 0),]

      if(!is.null(mask)){
        X_2D <- points_2D[mask,]
        X_garbage_2D <- points_2D[!mask, ]
        X_labels <- X_labels[mask]
      }

      colors_points <- rainbow(length(table(X_labels)))
      colors_labels <- sapply(X_labels, function(x){colors_points[x]})

      plot(points_2D[,1],points_2D[,2], col=colors_labels, pch=as.character(X_labels),...)
      points(boundary_ids_2D[,1],boundary_ids_2D[,2], pch="x",col="black",cex=4)

      colors_path <- rainbow(length(unique(perturbed_path_labels)))
      colors_labels_ppath <- sapply(perturbed_path_labels, function(x){colors_path[abs(x)]})

      i <- 1
      for(el in (unique(perturbed_path_labels))){
        lines(ppath_2D[which(perturbed_path_labels==el),1], ppath_2D[which(perturbed_path_labels==el),2],lwd=3,col=colors_labels_ppath[which(perturbed_path_labels == el)],type="o",pch=i)
        i <- i + 1
      }
      if(!is.null(mask)){
        points(X_garbage_2D[,1],X_garbage_2D[,2], col="gray", pch="x")
      }
    }
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
#' @examples
#' load("~/spathial/inst/extdata/X.rda")
#' load("~/spathial/inst/extdata/X_labels.rda")
#' boundaryRes <- spathialBoundaryIds(X, X_labels, mode=2, from=3, to=6)
#' NC <- 20
#' spathialRes <- spathialWay(boundaryRes$X, boundaryRes$X_labels, boundaryRes$boundary_ids, NC, neighbors = 1)
#' statistics <- spathialStatistics(spathialRes)
#' @export
spathialStatistics <- function(spathial_res){

  #Mean of z scores (with Fisher transformation)
  if(!is.null(spathial_res$perturbed_paths)){
    z_scores_perturbed_paths <- lapply(spathial_res$perturbed_paths, function(x){
      lapply(x, function(y){
        z_scores <- apply(y, 2, function(z){
          if(sd(z) == 0){
            return(0)
          }else{
            DescTools::FisherZ(cor(z, c(1:length(z))))
          }
        })
      })
    })

    count <- 0
    sum <- array(data=0, dim=c(ncol(spathial_res$ppath)))
    for(i in (1:length(z_scores_perturbed_paths))){
      A <- z_scores_perturbed_paths[[i]]
      for(j in (1:length(A))){
        B <- A[[j]]
        for(k in (1:length(B))){
          sum[k] <- sum[k] +  B[k]
        }
        count <- count + 1
      }
    }
    correlations <- sapply(sum, function(x){
      z_avg <- x/count
      return(DescTools::FisherZInv(z_avg))
    })
    correlations <- as.list(correlations)
    names(correlations) <- colnames(spathial_res$ppath)
    fisher<-unlist(correlations)

    #Mean of ranks
    ranks <- lapply(spathial_res$perturbed_path, function(x){
      lapply(x, function(y){
        corr <- apply(y, 2, function(z){
          if(sd(z) == 0){
            return(0)
          }else{
            cor(z, c(1:length(z)))
          }
        })
        return(rank(-corr))
      })
    })

    count <- 0
    sum <- array(data=0, dim=c(ncol(spathial_res$ppath)))
    for(i in (1:length(ranks))){
      A <- ranks[[i]]
      for(j in (1:length(A))){
        B <- A[[j]]
        for(k in (1:length(B))){
          sum[k] <- sum[k] +  B[k]
        }
        count <- count + 1
      }
    }

    ranks <- sapply(sum, function(x){
      rank_avg <- x/count
    })
    ranks <- as.list(ranks)
    names(ranks) <- colnames(spathial_res$ppath)
    ranks<-unlist(ranks)

  }else{ # When a single path is present
    correlations <- sapply(spathial_res$ppath, function(x){
      if(sd(x) == 0){
        return(0)
      }else{
        cor(x, c(1:length(x)))
      }
    })
    correlations<-unlist(correlations)
    names(correlations)<-colnames(spathial_res$ppath)
    ranks <- rank(-correlations)
    names(ranks)<-colnames(spathial_res$ppath)
  }

  outlist<-list(
    correlations=correlations,
    ranks=ranks
  )
  return(outlist)
}

