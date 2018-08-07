#### Class!! ---------------------------------------------------------------

impute_class <- setClass("impute_class", slots = c(imputations = "data.frame",
                                                   number_of_observations,
                                                   mcall = "call",
                                                   number_of_imputations = "integer",
                                                   censoring_type = "character"))

setMethod(f = "print", signature = "impute_class",
          function(imputations, mcall, number_of_imputations, censoring_type, ...){
            #stopifnot(inherits(x, what = "imputed_class"))
            
            #z <- x
            
            cl <- noquote(capture.output(mcall@mcall)) #call
            m <- number_of_imputations@number_of_imputations
            cens_type <- censoring_type@censoring_type
            n <- number_of_observations@number_of_observations
            
            cat("\n", 
                sprintf("Call: \n %s\n\n", cl),
                sprintf("Number of observations: %s\n", n),
                paste("Imputed", m, cens_type, "censored values") )
            
            
          })
