#' Find the Optimal Rerandomization Design Under the Tail and Kurtosis Approximation
#' 
#' Finds the optimal rerandomization threshold based on a user-defined quantile
#' and kurtosis based on an approximation of tail standard errors
#' 
#' @param W_base_object			An object that contains the assignments to begin with sorted by imbalance.
#' @param estimator 			"linear" for the covariate-adjusted linear regression estimator (default).
#' @param q 					The tail criterion's quantile of MSE over z's. The default is 95\%. 
#' @param c_val					The c value used (see Equation 8 in the paper). The default is \code{NULL} corresponding to \code{qnorm(q)}.
#' @param skip_search_length	In the exhaustive search, how many designs are skipped? Default is 1 for 
#' 								full exhaustive search through all assignments provided for in \code{W_base_object}.
#' @param binary_search			If \code{TRUE}, a binary search is employed to find the optimal threshold instead of 
#' 								an exhaustive search. Default is \code{FALSE}.
#' @param excess_kurtosis_z		An estimate of the excess kurtosis in the measure on z. Default is 0.
#' @param use_frob_norm_sq_unbiased_estimator		If \code{TRUE}, this would use the debiased Frobenius norm estimator
#' 													instead of the naive. Default is \code{TRUE}.
#' @param frob_norm_sq_bias_correction_min_samples  The bias-corrected estimate suffers from high variance when there 
#' 													are not enough samples. Thus, we only implement
#' 													the correction beginning at this number of vectors. Default is 10 and
#' 													this parameter is only applicable if \code{use_frob_norm_sq_unbiased_estimator} is \code{TRUE}.
#' @param smoothing_degree  	The smoothing degree passed to \code{loess}.
#' @param smoothing_span    	The smoothing span passed to \code{loess}.
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
#'  design = optimal_rerandomization_tail_approx(W_base_obj, 
#' 				skip_search_length = 10)
#'  design
#' 	}
#' @export
optimal_rerandomization_tail_approx = function(
		W_base_object,
		estimator = "linear",
		q = 0.95,
		c_val = NULL,
		skip_search_length = 1,
		binary_search = FALSE,
		excess_kurtosis_z = 0,
		use_frob_norm_sq_unbiased_estimator = TRUE,
		frob_norm_sq_bias_correction_min_samples = 10,
		smoothing_degree = 1,
		smoothing_span = 0.1,
		dot_every_x_iters = 100){
	
	optimal_rerandomization_argument_checks(W_base_object, estimator, q)
	
	if (is.null(c_val)){
		c_val = qnorm(q)
	}
	
	n = W_base_object$n
	X = W_base_object$X
	W_base_sort = W_base_object$W_base_sort
	max_designs = W_base_object$max_designs
	imbalance_by_w_sorted = W_base_object$imbalance_by_w_sorted
	
	if (estimator == "linear"){
	  Xt = t(X)
	  XtXinv = solve(Xt %*% X)
	  XtXinv_eigen = eigen(XtXinv)
	  if (W_base_object$p == 1){
		  E = sqrt(XtXinv_eigen$values)
	  } else {
		  E = diag(sqrt(XtXinv_eigen$values))
	  }
	  
	  XtXinv_sqrt = XtXinv_eigen$vectors %*% E %*% solve(XtXinv_eigen$vectors)
	  X_orth = X %*% XtXinv_sqrt
	  X_orth_T = t(X_orth)
	  P = X %*% XtXinv %*% Xt
	  I = diag(n)
	  I_min_P = I - P
	}


	Q_primes = array(NA, max_designs)
	#diagram all terms
	balances = array(NA, max_designs)
	frob_norm_sqs = array(NA, max_designs)
	tr_gds = array(NA, max_designs)
	tr_d_sqs = array(NA, max_designs)
	r_i_sqs = array(NA, max_designs)
	
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
	    balances[s] = tr(X_orth_T %*% Sigma_W %*% X_orth)
		if (use_frob_norm_sq_unbiased_estimator){
			frob_norm_sqs[s] = #see Appendix 5.9	 
				frob_norm_sq_debiased(Sigma_W, i, n, frob_norm_sq_bias_correction_min_samples)	-
				frob_norm_sq_debiased_times_matrix(Sigma_W, I_min_P, i, n, frob_norm_sq_bias_correction_min_samples)
		} else {
			frob_norm_sqs[s] = frob_norm_sq(I_min_P %*% Sigma_W)
		}

	    tr_gds[s] = tr(G %*% D) / n
	    tr_d_sqs[s] = tr(D^2) / n^2
	    r_i_sqs[s] = sum((diag(G) + 2 * diag(D) / n)^2)
	    Q_primes[s] = balances[s] +
	      c_val * sqrt(
		      2 * frob_norm_sqs[s] + 
	          8 * tr_gds[s] +
	          8 * tr_d_sqs[s] + 
			  excess_kurtosis_z * r_i_sqs[s]
	      )
	  
	  }
	}
	cat("\n")
	
	s_star = NULL
	Q_star = Inf
	for (s in seq(from = 1, to = max_designs, by = skip_search_length)){
	  if (Q_primes[s] < Q_star){
	    Q_star = Q_primes[s]
	    s_star = s
	  }
	}
	
	#now do the smoothing
	smoothing_fit = loess(Q_primes ~ imbalance_by_w_sorted, degree = smoothing_degree, span = smoothing_span)
	Q_primes_smoothed = predict(smoothing_fit, data.frame(X = imbalance_by_w_sorted))
	
	s_star_smoothed = NULL
	Q_star_smoothed = Inf
	for (s in seq(from = 1, to = max_designs, by = skip_search_length)){
		if (Q_primes_smoothed[s] < Q_star_smoothed){
			Q_star_smoothed = Q_primes_smoothed[s]
			s_star_smoothed = s
		}
	}
	
	all_data_from_run = data.frame(
	  imbalance_by_w_sorted = imbalance_by_w_sorted, 
	  Q_primes = Q_primes,
	  Q_primes_smoothed = Q_primes_smoothed,
	  balances = balances,
	  frob_norm_sqs = frob_norm_sqs,
	  tr_gds = tr_gds,
	  tr_d_sqs = tr_d_sqs,
	  r_i_sqs = r_i_sqs
	)

	ll = list(
		type = "approx",
		use_frob_norm_sq_unbiased_estimator = use_frob_norm_sq_unbiased_estimator,
		frob_norm_sq_bias_correction_min_samples = frob_norm_sq_bias_correction_min_samples,
		smoothing_degree = smoothing_degree,
		smoothing_span = smoothing_span,
		estimator = estimator,
		q = q,
		W_base_object = W_base_object,
		W_star = W_base_sort[1 : s_star, ],
		W_star_size = s_star,
		a_star = imbalance_by_w_sorted[s_star],
		a_stars = imbalance_by_w_sorted[1 : s_star],
		W_star_size_smoothed = s_star_smoothed,
		a_star_smoothed = imbalance_by_w_sorted[s_star_smoothed],
		a_stars_smoothed = imbalance_by_w_sorted[1 : s_star_smoothed],
		all_data_from_run = all_data_from_run,
		Q_star = Q_star
	)
	class(ll) = "optimal_rerandomization_obj"
	ll
}

