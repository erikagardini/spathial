\name{spathial_way}
\alias{spathial_way}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
spathial_way(X, boundaries, NC, prefiltering)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{boundaries}{
%%     ~~Describe \code{boundaries} here~~
}
  \item{NC}{
%%     ~~Describe \code{NC} here~~
}
  \item{prefiltering}{
%%     ~~Describe \code{prefiltering} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, boundaries, NC, prefiltering) 
{
    if (prefiltering) {
        prefiltered <- rkm_prefilter(X, boundary_ids, plot_ax = TRUE)
        X <- prefiltered$X_filtered
        boundary_ids <- prefiltered$boundary_ids_filtered
        X_g <- prefiltered$X_garbage
        rm(prefiltered)
    }
    waypoint_ids <- initMedoids(X, NC, "kpp", boundary_ids)
    waypoint_ids <- c(boundary_ids[1], waypoint_ids, boundary_ids[2])
    init_W <- X[waypoint_ids, ]
    s_span <- pracma::logspace(5, -5, n = NC)
    s_span <- c(s_span, 0)
    models <- list()
    for (i in 1:length(s_span)) {
        s <- s_span[i]
        W <- rkm(X, init_W, s, plot_ax = TRUE)
        init_W <- W
        models[[as.character(s)]] <- W
    }
    W_dst_var <- rkm_MS_pathvar(models, s_span, X)
    s_elb_id <- find_elbow(cbind(s_span, W_dst_var))
    return(models[[s_elb_id]])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
