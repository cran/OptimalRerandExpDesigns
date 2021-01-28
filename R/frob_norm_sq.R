#' Naive Frobenius Norm Squared
#' 
#' Compute naive / vanilla squared Frobenius Norm of matrix A
#' 
#' @param A			The matrix of interest
#' @return 			The Frobenius Norm of \code{A} squared.
#'
#' @author Adam Kapelner
#' @export
frob_norm_sq = function(A){sum(A^2)}

#' Debiased Frobenius Norm Squared Var-Cov matrix
#' 
#' Compute debiased Frobenius Norm of matrix Sigmahat (Appendix 5.8). Note
#' that for S <= 2, it returns the naive estimate.
#' 
#' @param Sigmahat  The var-cov matrix of interest
#' @param s         The number of vectors \code{Sigmahat} was generated from
#' @param n         The length of each vector
#' @param frob_norm_sq_bias_correction_min_samples  This estimate suffers from high variance when there are not enough samples. Thus, we only implement
#' the correction beginning at this number of samples otherwise we return the naive estimate. Default is 10.
#' @return 			The unbiased estimate of the Frobenius Norm of a variance-covariance matrix squared.
#'
#' @author Adam Kapelner
#' @export
frob_norm_sq_debiased = function(Sigmahat, s, n, frob_norm_sq_bias_correction_min_samples = 10){
  frob_norm_sq_naive = frob_norm_sq(Sigmahat)
  if (s <= frob_norm_sq_bias_correction_min_samples){
    frob_norm_sq_naive
  } else {
    s / (s - 1) * frob_norm_sq_naive - n^2 / (s - 1)
  }
}

#' Debiased Frobenius Norm Squared Constant Times Var-Cov matrix
#' 
#' Compute debiased Frobenius Norm of matrix P times Sigmahat (Appendix 5.9). Note
#' that for S <= 2, it returns the naive estimate.
#' 
#' @param Sigmahat  The var-cov matrix of interest
#' @param A         The matrix that multiplies Sigmahat
#' @param s         The number of vectors \code{Sigmahat} was generated from
#' @param n         The length of each vector
#' @param frob_norm_sq_bias_correction_min_samples  This estimate suffers from high variance when there are not enough samples. Thus, we only implement
#' the correction beginning at this number of samples otherwise we return the naive estimate. Default is 10.
#' @return 			The unbiased estimate of the Frobenius Norm of \code{A} times a variance-covariance matrix quantity squared.
#'
#' @author Adam Kapelner
#' @export
frob_norm_sq_debiased_times_matrix = function(Sigmahat, A, s, n, frob_norm_sq_bias_correction_min_samples = 10){
  frob_norm_sq_naive = frob_norm_sq(A %*% Sigmahat)
  if (s <= frob_norm_sq_bias_correction_min_samples){
    frob_norm_sq_naive
  } else {
    from_norm_sq_A_Sigmahat = 0
    for (i in 1 : n){
      from_norm_sq_A_Sigmahat = from_norm_sq_A_Sigmahat + t(A[i, ]) %*% Sigmahat %*% A[i, ]
    } 
    s / (s - 1) * frob_norm_sq_naive - n / (s - 1) * from_norm_sq_A_Sigmahat
  }
}