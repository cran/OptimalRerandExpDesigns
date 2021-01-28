#' A tool to find the optimal rerandomization threshold in non-sequential experiments
#'
#' @name 		OptimalRerandExpDesigns
#' @docType 	package
#' @title 		Optimal Rerandomization Threshold Search for Experimental Design
#' @author 		Adam Kapelner \email{kapelner@@qc.cuny.edu}
#' @references 	Kapelner, A
#' @importFrom  stats loess na.omit predict qchisq qnorm quantile rbinom sd uniroot var
#' @importFrom  ggplot2 aes geom_histogram geom_line geom_vline ggplot ggtitle scale_x_log10 scale_y_log10
#' @importFrom  GreedyExperimentalDesign initGreedyExperimentalDesignObject resultsGreedySearch compute_objective_val
#' @importFrom  momentchi2 hbe
##### Run "library(roxygen2); roxygenise("OptimalRerandExpDesigns", clean = TRUE)" to regenerate all Rd files and NAMESPACE and DESCRIPTION file
##### but make sure you are in the root directory of the project
NULL