% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimize_rerand_gaussian.R
\name{optimal_rerandomization_normality_assumed}
\alias{optimal_rerandomization_normality_assumed}
\title{Find the Optimal Rerandomization Design Under the Gaussian Approximation}
\usage{
optimal_rerandomization_normality_assumed(
  W_base_object,
  estimator = "linear",
  q = 0.95,
  skip_search_length = 1,
  dot_every_x_iters = 100
)
}
\arguments{
\item{W_base_object}{An object that contains the assignments to begin with sorted by}

\item{estimator}{"linear" for the covariate-adjusted linear regression estimator (default).}

\item{q}{The tail criterion's quantile of MSE over z's. The default is 95\%.}

\item{skip_search_length}{In the exhaustive search, how many designs are skipped? Default is 1 for 
full exhaustive search through all assignments provided for in \code{W_base_object}.}

\item{dot_every_x_iters}{Print out a dot every this many iterations. The default is 100. Set to
\code{NULL} for no printout.}
}
\value{
A list containing the optimal design threshold, strategy, and
								other information.
}
\description{
Finds the optimal rerandomization threshold based on a user-defined quantile
and a function that generates the non-linear component of the response
}
\examples{
 \donttest{
 n = 100
 p = 10
 X = matrix(rnorm(n * p), nrow = n, ncol = p)
 X = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})
 S = 25000
 
 W_base_obj = generate_W_base_and_sort(X, max_designs = S)
 design = optimal_rerandomization_normality_assumed(W_base_obj, 
				skip_search_length = 10)
 design
	}
}
\author{
Adam Kapelner
}
