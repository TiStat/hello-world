# rm(list = ls())
# setwd("C:/Users/Petro Ck/Google Drive/R-Project/R")



# Create summary method for "impute_class" class ------------------------------------------------
summary.impute_class <- function(x) {
  #stopifnot(inherits(x = result, what = "imputed_class"))
  
  # Referencing to our object
  # z <- object
  
  m <- x$number_of_imputations #number of imputations
  n <- x$number_of_observations
 #  v <- z$variance #discrepancy
 #  n <- z$number_of_observations
 #  aic <- z$AIC
  mcall <- noquote(capture.output(x$mcall)) #call
  cens_type <- x$censoring_type
  
  cat("\n", 
            sprintf("Call: \n %s\n\n", mcall),
            sprintf("Number of observations: %s\n", n), 
            sprintf("Number of imputations: %s", m))
                
                
  
}

# TO BE WRITTEN IN MAIN FUNCTION!! ---------------------------------------------------------------

# All of the above

print.impute_class <- function(x, ...){
 # stopifnot(inherits(x, what = "imputed_class"))
  
  #z <- result
  
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



### Some testing -----------------------------------------------
print.impute_class <- function(object){
  cl <- oldClass(x)
  oldClass(x) <- cl[cl != "impute_class"]
  NextMethod("print")
  invisible(x)
  # stopifnot(inherits(object, what = "imputed_class"))
  # 
  # z <- object
  # 
  # mcall <- noquote(capture.output(z$mcall)) #call
  # m <- z$number_of_imputations
  # cens_type <- z$censoring_type
  # 
  # cat("\n", 
  #     sprintf("Call: \n %s\n\n", mcall),
  #     sprintf("Number of observations: %s\n", n),
  #     paste("Imputed", m, cens_type, "censored values") )
  
  
}




print.impute_class <- function(x, ...){
  cl <- oldClass(x)
  oldClass(x) <- cl[cl != "impute_class"]
  NextMethod("print")
  invisible(x)
}




### Example ------------------------------------------------------------
  myMatrix <- function(x) {
    if(!is.matrix(x)){
      stop("x is not a matrix")
    }
    x <- x*2
    class(x) <- "myMatrix"  ## wichtig!!
    return(x)
  }

  print.myMatrix <- function(x) {
    n <- nrow(x)
    for(i in seq_len(n)) {
      cat(paste("This is row", i, "\t: " ))
      cat(x[i,], "\n")
    }
  }

a <- myMatrix(matrix(1:16, ncol=4))
a # wunderbar




