
#' Plots a summary of the imbalances in a \code{W_base_object} object
#' 
#' @param x			The \code{W_base_object} object to be summarized in the plot
#' @param ...		\code{title}, \code{subtitle}, \code{xlab}, \code{bins} can 
#' 					be specified here to be passed to the ggplot plotting function.
#' 					Also \code{log10} can be set to \code{FALSE} to not log the x-axis.
#' @return 			No return value, called for side effects
#' 
#' @author 			Adam Kapelner
#' @method plot W_base_object
#' @export
plot.W_base_object = function(x, ...){
  dots = list(...)
  if (is.null(dots$title)){
    title = "Density of Imbalances in Base Strategy"
  } else {
    title = dots$title
  }
  if (is.null(dots$subtitle)){
    subtitle = ""
  } else {
    subtitle = dots$subtitle
  }
  if (is.null(dots$xlab)){
    xlab = x$imbalance_function
  } else {
    xlab = dots$xlab
  }
  if (is.null(dots$bins)){
    bins = x$max_designs / 10
  } else {
    bins = dots$bins
  }
  
  ggplot_obj = ggplot(data.frame(b = x$imbalance_by_w_sorted)) + 
    ggtitle(title, subtitle = subtitle) +
    xlab(xlab) +
    geom_histogram(aes(x = b), bins = bins)
  if (!isFALSE(dots$log10)){
    ggplot_obj = ggplot_obj + scale_x_log10()
  }
  
  plot(ggplot_obj)
}

#' Prints a summary of a \code{W_base_object} object
#' 
#' @param x			The \code{W_base_object} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' @return 			No return value, called for side effects
#' 
#' @author 			Adam Kapelner
#' @method print W_base_object
#' @export
print.W_base_object = function(x, ...){	
  cat("W base strategy with", x$max_designs, "assignments whose imbalances range from",
      round(min(x$imbalance_by_w_sorted), 3), "to", round(max(x$imbalance_by_w_sorted), 3), 
      "in", x$imbalance_function, "\n")
}

#' Prints a summary of a \code{W_base_object} object
#' 
#' @param object		The \code{W_base_object} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary W_base_object
#' @export
summary.W_base_object = function(object, ...){
  print(object, ...)
}


#' Prints a summary of a \code{optimal_rerandomization_obj} object
#' 
#' @param x			The \code{optimal_rerandomization_obj} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' @return 			No return value, called for side effects
#' 
#' @author 			Adam Kapelner
#' @method print optimal_rerandomization_obj
#' @export
print.optimal_rerandomization_obj = function(x, ...){	
  if (x$type == "exact") {
    cat("Optimal rerandomization found with", x$W_star_size_smoothed, "assignments whose imbalances are smaller\nthan",
        round(x$a_star_smoothed, 3), "in", x$imbalance_function, "using algorithm type", x$type, "(smoothed) at q =", x$q, "\n")} else {
    cat("Optimal rerandomization found with", x$W_star_size, "assignments whose imbalances are smaller\nthan",
        round(x$a_star, 3), "in", x$imbalance_function, "using algorithm type", x$type, "at q =", x$q, "\n")
  }
}

#' Prints a summary of a \code{optimal_rerandomization_obj} object
#' 
#' @param object		The \code{optimal_rerandomization_obj} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary optimal_rerandomization_obj
#' @export
summary.optimal_rerandomization_obj = function(object, ...){
  print(object, ...)
}

#' Plots a summary of a \code{optimal_rerandomization_obj} object
#' 
#' @param x			The \code{optimal_rerandomization_obj} object to be summarized in the plot
#' @param ...		The option \code{advanced = TRUE} can be passed here for optimal rerandomization 
#' 					results from algorithm type "approx" to see how all the terms in the criterion behave.
#' 					You can pass \code{s_min} which controls the minimum number of vectors the plot begins at. 
#' 					Below a certain number, the criterion is unstable.
#' 					Also, \code{title}, \code{subtitle}, \code{xlab} and \code{ylab} can be passed here.
#' @return 			No return value, called for side effects
#' 
#' @author 			Adam Kapelner
#' @method plot optimal_rerandomization_obj
#' @export
plot.optimal_rerandomization_obj = function(x, ...){
  dots = list(...)
  if (is.null(dots$title)){
    title = "Optimal rerandomization by Tail Criterion and All Terms"
  } else {
    title = dots$title
  }
  if (is.null(dots$subtitle)){
    subtitle = 
    if (x$type == "normal"){
      subtitle = paste("optimal indicated by green line:", x$W_star_size, "of", x$W_base_object$max_designs, "vectors")
    } else {
      subtitle = paste("optimal / smoothed optimal indicated by green / purple line:", x$W_star_size, "/", x$W_star_size_smoothed, "of", x$W_base_object$max_designs, "vectors")
    }
  } else {
    subtitle = dots$subtitle
  }
  if (is.null(dots$xlab)){
    xlab = x$imbalance_function
  } else {
    xlab = dots$xlab
  }
  if (is.null(dots$ylab)){
    ylab = paste("Relative MSE Tail at q =", x$q)
  } else {
    ylab = dots$ylab
  }
  if (is.null(dots$s_min)){
    s_min = 1
  } else {
    s_min = dots$s_min
  }
  max_designs = nrow(x$all_data_from_run)
  
  if (x$type == "approx" && isTRUE(dots$advanced)){
    plot(ggplot(na.omit(x$all_data_from_run[s_min : max_designs, ])) +
           ggtitle(title, subtitle = subtitle) +
           xlab(xlab) +
           ylab(ylab) +
           scale_x_log10() +
           scale_y_log10() +
           geom_line(aes(x = imbalance_by_w_sorted, y = Q_primes)) +
		   geom_line(aes(x = imbalance_by_w_sorted, y = Q_primes_smoothed), col = "purple") +
           geom_line(aes(x = imbalance_by_w_sorted, y = imbalance_by_w_sorted), col = "blue") +
           geom_line(aes(x = imbalance_by_w_sorted, y = frob_norm_sqs), col = "red") +
           geom_line(aes(x = imbalance_by_w_sorted, y = tr_gds), col = "orange") +
           geom_line(aes(x = imbalance_by_w_sorted, y = tr_d_sqs), col = "grey") + 
           geom_line(aes(x = imbalance_by_w_sorted, y = r_i_sqs), col = "yellow") +			
           geom_vline(xintercept = log(x$a_star), col = "green"))
  } else if (x$type == "exact" || x$type == "approx") {
    plot(ggplot(na.omit(x$all_data_from_run[s_min : max_designs, ])) +
           ggtitle(title, subtitle = subtitle) +
           xlab(xlab) +
           ylab(ylab) +
           scale_x_log10() +
           scale_y_log10() +
           geom_line(aes(x = imbalance_by_w_sorted, y = Q_primes)) +
		   geom_line(aes(x = imbalance_by_w_sorted, y = Q_primes_smoothed), col = "purple") +
           geom_vline(xintercept = x$a_star, col = "green") +
           geom_vline(xintercept = x$a_star_smoothed, col = "purple"))
  } else {
    plot(ggplot(na.omit(x$all_data_from_run[s_min : max_designs, ])) +
           ggtitle(title, subtitle = subtitle) +
           xlab(xlab) +
           ylab(ylab) +
           scale_x_log10() +
           scale_y_log10() +
           geom_line(aes(x = imbalance_by_w_sorted, y = Q_primes)) +
           geom_vline(xintercept = x$a_star, col = "green"))
  }
  
}