#' z2p
#'
#' This function gives a gaussian p-value corresponding to the provided Z-score
#'
#' @param z a Z score
#' @return a p-value
#' @examples
#' z<-1.96
#' z2p(z)
#' @export
z2p <- function(z) {
    pnorm(abs(z), lower.tail = FALSE) * 2
}

#' p2z
#'
#' This function gives a gaussian Z-score corresponding to the provided p-value
#' Careful: sign is not provided
#'
#' @param p a p-value
#' @return z a Z score
#' @examples
#' p<-0.05
#' p2z(p)
#' @export
p2z <- function(p) {
    qnorm(p/2, lower.tail = FALSE)
}

#' Stouffer integration of Z scores
#'
#' This function gives a gaussian Z-score corresponding to the provided p-value
#' Careful: sign is not provided
#'
#' @param x a vector of Z scores
#' @return Z an integrated Z score
#' @examples
#' zs<-c(1,3,5,2,3)
#' stouffer(zs)
#' @export
stouffer <- function(x) {
    Z <- sum(x)/sqrt(length(x))
    return(Z)
}

#' Weighted Stouffer integration of Z scores
#'
#' This function gives a gaussian Z-score corresponding to the provided p-value
#' Careful: sign is not provided
#'
#' @param x a vector of Z scores
#' @param w weight for each Z score
#' @return Z an integrated Z score
#' @examples
#' zs<-c(1,-3,5,2,3)
#' ws<-c(1,10,1,2,1)
#' wstouffer(zs,ws)
#' @export
wstouffer <- function(x, w) {
    Z <- sum(x * w)/sqrt(sum(w^2))
    return(Z)
}


#' Fisher integration of p-values
#'
#' This function applies the Fisher integration of pvalues
#'
#' @param ps a vector of p-values
#' @return p.val an integrated p-value
#' @examples
#' ps<-c(0.01,0.05,0.03,0.2)
#' fisherp(ps)
#' @export
fisherp <- function(ps) {
    Xsq <- -2 * sum(log(ps))
    p.val <- pchisq(Xsq, df = 2 * length(ps),
                    lower.tail = FALSE)
    # p<-c(Xsq = Xsq, p.value = p.val)
    return(p.val)
}


#' Slice
#'
#' This function prints a slice of a matrix
#'
#' @param matrix A matrix
#' @return prints it
#' @examples
#' set.seed(1)
#' example<-matrix(rnorm(1000),nrow=100,ncol=10)
#' slice(example)
#' @export
slice <- function(matrix) {
    if (nrow(matrix) < 5) {
        stop("Input matrix has less than 5 rows")
    }
    if (ncol(matrix) < 5) {
        stop("Input matrix has less than 5 columns")
    }
    print(matrix[1:5, 1:5])
}



#' Convert correlation coefficient to p-value
#'
#' This functions converts an R value from a correlation calculation into a
#' p-value, using a T distribution with the provided number of samples N minus
#' 2 degrees of freedom
#' @param r a correlation coefficient
#' @param N a number of samples
#' @return p a p-value
#' @examples
#' set.seed(1)
#' a<-rnorm(1000)
#' b<-a+rnorm(1000,sd=10)
#' r<-cor(a,b,method='pearson')
#' corr2p(r,N=length(a))
#' @export
corr2p <- function(r, N) {
    # Get the t-distribution value
    t <- r/sqrt((1 - r^2)/(N - 2))
    # t follows the t distribution with
    # df=N?2
    p <- 1 - pt(abs(t), N - 2)
    return(p)
}

#' Convert p-value to correlation coefficient
#'
#' This functions converts an p-value into a the corresponding correlation
#' coefficient, using a T distribution with the provided number of samples N
#' minus 2 degrees of freedom
#' @param p a p-value
#' @param N a number of samples
#' @return r a correlation coefficient
#' @examples
#' N<-100
#' p<-0.05
#' p2corr(p,N)
#' @export
p2corr <- function(p, N) {
    # Get the t-value
    t <- abs(qt(p, N - 2))

    # Convert to correlation
    r <- sqrt(t^2/(N - 2 + t^2))

    return(r)
}

#' kmgformat - Nice Formatting of Numbers
#'
#' This function will convert thousand numbers to K, millions to M, billions
#' to G, trillions to T, quadrillions to P
#'
#' @param input A vector of values
#' @param roundParam How many decimal digits you want
#' @return A character vector of formatted numebr names
#' @examples
#' # Thousands
#' set.seed(1)
#' a<-runif(1000,0,1e4)
#' plot(a,yaxt='n')
#' kmg<-kmgformat(pretty(a))
#' axis(2,at=pretty(a),labels=kmg)
#'
#' # Millions to Billions
#' set.seed(1)
#' a<-runif(1000,0,1e9)
#' plot(a,yaxt='n',pch=20,col="black")
#' kmg<-kmgformat(pretty(a))
#' axis(2,at=pretty(a),labels=kmg)
#' @export
kmgformat <- function(input, roundParam = 1) {
    signs <- sign(input)
    signs[signs == 1] <- ""
    signs[signs == -1] <- "-"
    absinput <- abs(input)
    output <- c()
    for (i in absinput) {
        if (i < 1000) {
            output <- c(output, i)
        } else if (i < 1e+06) {
            i <- round(i/1000, roundParam)
            i <- paste0(i, "K")
            output <- c(output, i)
        } else if (i < 1e+09) {
            i <- round(i/1e+06, roundParam)
            i <- paste0(i, "M")
            output <- c(output, i)
        } else if (i < 1e+12) {
            i <- round(i/1e+09, roundParam)
            i <- paste0(i, "G")
            output <- c(output, i)
        } else if (i < 1e+15) {
            i <- round(i/1e+12, roundParam)
            i <- paste0(i, "T")
            output <- c(output, i)
        } else if (i < 1e+18) {
            i <- round(i/1e+15, roundParam)
            i <- paste0(i, "P")
            output <- c(output, i)
        } else {
            output <- c(output, i)
        }
    }
    output <- paste0(signs, output)
    return(output)
}

