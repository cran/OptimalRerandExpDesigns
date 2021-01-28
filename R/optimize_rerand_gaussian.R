
#' Find the Optimal Rerandomization Design Under the Gaussian Approximation
#' 
#' Finds the optimal rerandomization threshold based on a user-defined quantile
#' and a function that generates the non-linear component of the response
#' 
#' @param W_base_object			An object that contains the assignments to begin with sorted by 
#' @param estimator 			"linear" for the covariate-adjusted linear regression estimator (default).
#' @param q 					The tail criterion's quantile of MSE over z's. The default is 95\%. 
#' @param skip_search_length	In the exhaustive search, how many designs are skipped? Default is 1 for 
#' 								full exhaustive search through all assignments provided for in \code{W_base_object}.
#' @param dot_every_x_iters		Print out a dot every this many iterations. The default is 100. Set to
#' 								\code{NULL} for no printout.
#' @return 						A list containing the optimal design threshold, strategy, and
#' 								other information.
#' 
#' @author Adam Kapelner
#' @examples
#'  \donttest{
#'  n = 100
#'  p = 10
#'  X = matrix(rnorm(n * p), nrow = n, ncol = p)
#'  X = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})
#'  S = 25000
#'  
#'  W_base_obj = generate_W_base_and_sort(X, max_designs = S)
#'  design = optimal_rerandomization_normality_assumed(W_base_obj, 
#' 				skip_search_length = 10)
#'  design
#' 	}
#' @export
optimal_rerandomization_normality_assumed = function(
  W_base_object,
  estimator = "linear",
  q = 0.95,
  skip_search_length = 1,
  dot_every_x_iters = 100){
  optimal_rerandomization_argument_checks(W_base_object, estimator, q)
  
  n = W_base_object$n
  X = W_base_object$X
  W_base_sort = W_base_object$W_base_sort
  max_designs = W_base_object$max_designs
  imbalance_by_w_sorted = W_base_object$imbalance_by_w_sorted
  
  if (estimator == "linear"){
    Xt = t(X)
    XtXinv = solve(Xt %*% X)
    P = X %*% XtXinv %*% Xt
    I = diag(n)
    I_min_P = I - P
  }
  
  s_star = NULL
  Q_star = Inf
  Q_primes = array(NA, max_designs)
  
  w_w_T_running_sum = matrix(0, n, n)
  if (estimator == "linear"){
    w_w_T_P_w_w_T_running_sum = matrix(0, n, n)
  }
  ss = seq(from = 1, to = max_designs, by = skip_search_length)
  for (i in 1 : length(ss)){
	s = ss[i]
    if (!is.null(dot_every_x_iters)){
      if (i %% dot_every_x_iters == 0){
        cat(".")
      }
    }
    w_s = W_base_sort[s, , drop = FALSE]
    w_s_w_s_T = t(w_s) %*% w_s
    w_w_T_running_sum = w_w_T_running_sum + w_s_w_s_T
    Sigma_W = 1 / i * w_w_T_running_sum
    if (estimator == "linear"){			
      w_w_T_P_w_w_T_running_sum = w_w_T_P_w_w_T_running_sum + w_s_w_s_T %*% P %*% w_s_w_s_T
      D = 1 / i * w_w_T_P_w_w_T_running_sum
      G = I_min_P %*% Sigma_W %*% I_min_P
      eigenvalues = eigen(G + 2 / n * D)$values
	  
	  #sometimes the hall_buckley_eagleson method fails. I think it's due to numerical instability in the algorithm. 
	  #We will simply skip these.
	  tryCatch({Q_primes[s] = hall_buckley_eagleson_inverse_cdf(eigenvalues, q, n)}, error = function(e){})	
    }
	if (!is.na(Q_primes[s])){
		if (Q_primes[s] < Q_star){
			Q_star = Q_primes[s]
			s_star = s
		}		
	}

  }
  cat("\n")
  
  all_data_from_run = data.frame(
    imbalance_by_w_sorted = imbalance_by_w_sorted, 
    Q_primes = Q_primes
  )
  
  ll = list(
    type = "normal",
    estimator = estimator,
    q = q,
    W_base_object = W_base_object,
    W_star = W_base_sort[1 : s_star, ],
    W_star_size = s_star,
    a_star = imbalance_by_w_sorted[s_star],
    a_stars = imbalance_by_w_sorted[1 : s_star],
    all_data_from_run = all_data_from_run,
    Q_star = Q_star
  )
  class(ll) = "optimal_rerandomization_obj"
  ll
}