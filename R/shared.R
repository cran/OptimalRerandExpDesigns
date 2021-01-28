#optimal_rerandomization_sets = function(X, 
#                                        max_designs = 25000,
#                                        objective = "mahal_dist",
#                                        estimator = "linear",
#                                        c_grid = seq(from = 1, to = 5, by = 0.5),
#                                        kappa_z_grid = seq(from = -1, to = 1, by = 0.5),
#                                        ...){
#  
#  if (estimator == "linear"){
#    all_results_a_star = array(NA, c(length(c_grid), length(kappa_z_grid)))
#    dimnames(all_results_a_star) = list(c_grid, kappa_z_grid)
#    all_results_Q_star = array(NA, c(length(c_grid), length(kappa_z_grid)))
#    dimnames(all_results_Q_star) = list(c_grid, kappa_z_grid)
#    all_results_W_star_size = array(NA, c(length(c_grid), length(kappa_z_grid)))
#    dimnames(all_results_W_star_size) = list(c_grid, kappa_z_grid)
#    for (i_c in 1 : length(c_grid)){
#      c_val = c_grid[i_c]
#      for (i_k in 1 : length(kappa_z_grid)){
#        kappa_z = kappa_z_grid[i_k]
#        
#        results = optimal_rerandomization(
#          X = X,
#          max_designs = max_designs,
#          objective = objective,
#          estimator = estimator,
#          c_val = c_val,
#          kappa_z = kappa_z
#        )
#        
#        all_results_a_star[i_c, i_k] = results$a_star
#        all_results_Q_star[i_c, i_k] = results$Q_star
#        all_results_W_star_size[i_c, i_k] = results$W_star_size
#      }
#    }
#  }
#  list(all_results_a_star = all_results_a_star, all_results_Q_star = all_results_Q_star, all_results_W_star_size = all_results_W_star_size)
#  load( "different_c_kappa.RData")
#  #save(ll, file = "different_c_kappa.RData")
#}





#' Generate Base Assignments and Sorts
#' 
#' Generates the base vectors to be used when locating the optimal rerandomization threshold
#' 
#' @param X						The data as an \eqn{n \times p} matrix. 
#' @param max_designs 			The maximum number of designs. Default is 25,000.
#' @param imbalance_function 	A string indicating the imbalance function. Currently,
#' 								"abs_sum_difference" and "mahal_dist" are the options with the
#' 								latter being the default.
#' @param r 					An experimental feature that adds lower imbalance vectors
#' 								to the base set using the \code{GreedyExperimentalDesign}
#' 								package. This controls the number of vectors to search through
#' 								on each iteration.
#' @param max_max_iters			An experimental feature that adds lower imbalance vectors
#' 								to the base set using the \code{GreedyExperimentalDesign}
#' 								package. The maximum number of iterations to use for the greedy search.
#' @return 						A list including all arguments plus a matrix \code{W_base_sorted}
#' 								whose \code{max_designs} rows are \code{n}-length allocation vectors
#' 								and the allocation vectors are in 
#'
#' 
#' @author Adam Kapelner
#' @examples
#'  n = 100
#'  p = 10
#'  X = matrix(rnorm(n * p), nrow = n, ncol = p)
#'  X = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})
#'  S = 1000
#'  
#'  W_base_obj = generate_W_base_and_sort(X, max_designs = S)
#'  W_base_obj
#' @export
generate_W_base_and_sort = function(X,
                                    max_designs = 25000,
                                    imbalance_function = "mahal_dist",
                                    r = 0,
                                    max_max_iters = 5
){
  n = nrow(X)
  p = ncol(X)
  #rewrite it as standardized
  X = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})
  
  W_base = complete_randomization_with_forced_balance_plus_one_min_one(n, max_designs)  
  
  if (r > 0){
    
    rd = initGreedyExperimentalDesignObject(X, r, wait = TRUE, objective = "mahal_dist", num_cores = 4)
    res = resultsGreedySearch(rd, max_vectors = r)
    all_vecs = t(res$ending_indicTs)
    all_vecs[all_vecs == 0] = -1
    W_base = rbind(all_vecs, W_base)
    
    for (max_iters in 1 : max_max_iters){
      rd = initGreedyExperimentalDesignObject(X, r, wait = TRUE, objective = "mahal_dist", num_cores = 4, max_iters = max_iters)
      res = resultsGreedySearch(rd, max_vectors = r)
      all_vecs = t(res$ending_indicTs)
      all_vecs[all_vecs == 0] = -1
      W_base = rbind(all_vecs, W_base)
    }
    
	all_mahal_obj_sort_log = log(sort(apply(W_base, 1, function(w){compute_objective_val(X, w, objective = "mahal_dist")})))
    
    
    ##now create a representative sample
    bin_width_in_log = 0.1
    tab = table(round(all_mahal_obj_sort_log / bin_width_in_log))
    
    num_per_bin_star = floor(uniroot(function(num_per_bin){
      sum(pmin(tab, num_per_bin)) - max_designs
    }, interval = c(1, max_designs))$root)
    sum(pmin(tab, num_per_bin_star))
    
    all_vecs_by_bin = list()
    for (i in 1 : nrow(W_base)){
      idx = as.character(round(all_mahal_obj_sort_log[i] / bin_width_in_log))
      if (is.null(all_vecs_by_bin[[idx]])){
        all_vecs_by_bin[[idx]] = matrix(NA, nrow = 0, ncol = n)
      }
      if (nrow(all_vecs_by_bin[[idx]]) < num_per_bin_star){
        all_vecs_by_bin[[idx]] = rbind(all_vecs_by_bin[[idx]], W_base[i, ])
      }
    }
    unlist(lapply(all_vecs_by_bin, nrow))
    W_base = df <- do.call("rbind", all_vecs_by_bin)
  }
  
  if (imbalance_function == "mahal_dist"){
    S_Xstd_inv = solve(var(X))
    imbalance_by_w = apply(W_base, 1, function(w){compute_objective_val_plus_one_min_one_enc(X, w, "mahal_dist", S_Xstd_inv)})
  } else {
    imbalance_by_w = apply(W_base, 1, function(w){compute_objective_val_plus_one_min_one_enc(X, w, imbalance_function)})
  }
  sorted_idx = sort(imbalance_by_w, index.return = TRUE)$ix
  
  ll = list(
    X = X,
    n = n,
	p = p,
    imbalance_function = imbalance_function,
    W_base_sorted = W_base[sorted_idx, ], 
    max_designs = nrow(W_base), 
    imbalance_by_w_sorted = imbalance_by_w[sorted_idx]
  )
  class(ll) = "W_base_object"
  ll
}



#' Returns the objective value given a design vector as well an an objective function.
#' This is code duplication since this is implemented within Java. This is only to be
#' run if...
#' 
#' @param X 		 	The n x p design matrix
#' @param indic_T		The n-length binary allocation vector
#' @param objective		The objective function to use. Default is \code{abs_sum_diff}.
#' @param inv_cov_X		Optional: the inverse sample variance covariance matrix. Use this
#' 						argument if you will be doing many calculations since passing this
#' 						in will cache this data.
#' @return 				A vector of computed objective values.
#' 
#' @author Adam Kapelner
#' @export
compute_objective_val_plus_one_min_one_enc = function(X, indic_T, objective = "abs_sum_diff", inv_cov_X = NULL){
  X_T = X[indic_T == 1, , drop = FALSE] #coerce as matrix in order to matrix multiply later
  X_C = X[indic_T == -1, , drop = FALSE] #coerce as matrix in order to matrix multiply later
  X_T_bar = colMeans(X_T)
  X_C_bar = colMeans(X_C)	
  
  if (objective == "abs_sum_diff"){
    s_j_s = apply(X, 2, sd)
    sum(abs((X_T_bar - X_C_bar) / s_j_s))
  } else if (objective == "mahal_dist"){
    #saves computation to pass it in if you're doing a lot of them in a row
    if (is.null(inv_cov_X)){
      inv_cov_X = solve(var(X))
    }	
    X_T_bar_minus_X_C_bar = as.matrix(X_T_bar - X_C_bar) #need to matricize for next computation
    as.numeric(t(X_T_bar_minus_X_C_bar) %*% inv_cov_X %*% X_T_bar_minus_X_C_bar)
  }
}

#compute trace of matrix
tr = function(A){sum(diag(A))}


#compute inverse CDF as a function of desired quantile using the hall-buckley-eagleson method
hall_buckley_eagleson_inverse_cdf = function(eigenvalues, q, sample_size, tol = 0.001){
  eigenvalues[eigenvalues <= 0] = .Machine$double.eps #only use non-zero eigenvalues
  if (any(eigenvalues > sample_size)){
    eigenvalues = c(sample_size) #theoretical upper limit
  }
  fun = function(x, eigenvalues){
    #cat("x =", x, "eigenvalues =", eigenvalues)
    hbe(coeff = eigenvalues, x = x) - q
  }
  uniroot(fun, eigenvalues, interval = c(tol, sample_size * qchisq(q, 1)) * 1.1, tol = tol)$root #1.1 is a fudge for numerical error.
}


#checks for illegal arguments
optimal_rerandomization_argument_checks = function(W_base_object, estimator, q){
  if (class(W_base_object) != "W_base_object"){
    stop("W_base_object must be class type \"W_base_object\".")
  }
  if (!(estimator %in% c("linear"))){
    stop("estimator must be \"linear\".")
  }
  if (q <= 0 || q >= 1){
    stop("quantile must be between 0 and 1.")
  }
}

#stupidity for CRAN!!!
utils::globalVariables(c(
	"all_mahal_obj_sort_log", 
	"b", 
	"imbalance_by_w_sorted", 
	"Q_primes", 
	"Q_primes_smoothed",
	"frob_norm_sqs",
	"tr_gds",
	"tr_d_sqs",
	"r_i_sqs"
))

		

