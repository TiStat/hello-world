# Create prin method for "imputed" class
#' Title
#'
#' @param x 
#'
#' @return What should be printed for "imputed" class object
#' @export
print.imputed <- function(x){
  
  mcall <- noquote(capture.output(x$mcall)) #call
  m <- x$number_of_imputations
  cens_type <- x$censoring_type
  n <- x$number_of_observations
  
  cat("\n", 
      cat("Call:"),
      mcall,
      sprintf("\n\n Number of observations: %s\n", n),
      paste("Imputed", m, cens_type, "censored values") )
}


# Create summary method for "imputed" class ------------------------------------------------
#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples 
summary.imputed <- function(x) {
  
  m <- x$number_of_imputations #number of imputations
  n <- x$number_of_observations
  #  v <- x$variance #discrepancy
  #  aic <- x$AIC
  mcall <- noquote(capture.output(x$mcall)) #call
  cens_type <- x$censoring_type
  
  cat("\n", 
      cat("Call:  \n", mcall),
      paste("\n Number of observations:", n),
      paste("\n Type of censoring:", cens_type),
      paste("\n Number of imputations:", m))
}