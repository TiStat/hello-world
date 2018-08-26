#' @title  Evaluate gamlss family functions from gamlss object
#' 
#' @description Function that, given a "gamlss" object, evaluates the
#'   distribution under predicted parameters of the provided dataframe
#'   "predictdata".
#'   CAUTION: Exactly ONE of the arguments x,q,p,n MUST be specfied!
#'   
#' @param object gamlss fit object
#' 
#' @param func character. "d", "p", "q", "r" for either density, distribution function,
#'   quantile or random data generation.
#'   
#' @param fitdata dataframe. Data used as input.
#' 
#' @param predictdata dataframe. Containing the observations for which the
#'   parameters are predicted.
#'   
#' @param x,q scalar numeric. Quantile value if density or probability function
#'   used respectively.
#'   
#' @param p scalar numeric. Probability value if quantile function used.
#' 
#' @param n scalar numeric. Number of observations if random generator function used.
#' 
#' @param ... argumenst to be passed to the called distributional function.
#'
#' @return Depending on the choice of func, the respective vector of
#'   (d)density-, (p)cdf.- , (q)quantile- or (r)random values is returned.
#'   
#' @export

family_fun <- function(object, func = c('d', 'p', 'q', 'r'), fitdata, predictdata ,p = NULL, q = NULL, x = NULL, n = NULL, ...) {
  
  func <- match.arg(func)
  
  if(!is.null(n)){
    if(func == 'r' & n%%nrow(predictdata) != 0){
      stop('Length of provided mu vector for rfamily is not a multiple of nrow(predictdata). 
            n > length(mu) implies, that the draws from multivariate 
            distributions are stacked. If n is not a multiple, the last vector that
            is to be stacked will be drawn from a shorter multivariate distribution.')
    }
  }
  
  # Find the correct family function to evaluate:
  fam_name <- object$family[1]
  f_fun <- paste(func, fam_name, sep= '')
  
  prd <- predictAll(
    object = object,
    type = "response",
    newdata = predictdata,
    data = fitdata
  )
  
  param <- list(
    mu = prd$mu,
    sigma = prd$sigma,
    nu = prd$nu,
    tau = prd$tau,
    x = x,
    q = q,
    p = p,
    n = n,
    ...
  )
  # Kick out NULL parameters:
  param <- param[!sapply(param, is.null)] 
  
  if (any(!names(param) %in% names(formals(f_fun)))) {
    stop("One of x,q,p,n,... arguments doesn't match with the distributional function 
         (e.g. dNO, pNO, qNO, rNO). See the family's documentation for admissable arguments.")
  }
  
  return(do.call(f_fun, param))
}

#' @title Inverse sampling - GAMLSS
#' 
#' @description Inverse sampling of censored variables to impute only valid
#'   observations, conditional on the respective fit.
#'   
#' @param object gamlss. Fitted model whose parameters are predicted for
#'   predictdata.
#'   
#' @param censtype character. Specifies the type of censoring (right/left/interval/missing).
#' 
#' @param predictdata dataframe. predict-data of the missing/censored
#'   observations, for which imputations are drawn.
#'   
#' @param fitdata data.frame. The orignal dataset upon which gamlss was fitted.
#' 
#' @param censor character. Name of the (damaged) column to be predicted on in
#'   predictdata. This is only required if censtype is NOT "missing".
#'   
#' @param intervalstart character. Name of the column of the interval's starting
#'   values. By convention, the starting duration in this column is assumed to
#'   be the time passed without failure, before entering the interval, in which
#'   the exact time of failure is unknown.
#'   
#' @param quantiles vector. Containing the quantiles to be evaluated in the
#'   conditoned distribution, i.e. conditoned on the parameters and the
#'   information contained in the censored value.
#'
#' @return Returns draws and quantiles.
#' 
#' @export

samplecensored <- function(object,
                           censtype,
                           predictdata,
                           fitdata,
                           censor,
                           quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
                           intervalstart = NULL) {
  
  # quantprob i.e. an auxilliary data frame with row-wise repeated quantiles vector:
  quantprob <- as.data.frame(matrix(rep(quantiles, times = nrow(predictdata)),
                                    byrow = TRUE, nrow = nrow(predictdata)))
  
  # WICKHAM STYLE: Make calls more readable! Always same arguments are passed
  f <- function(object, fitdata, predictdata) {
    g <- function(func, p = NULL, q = NULL, x = NULL, n = NULL) {
      family_fun(object, func, fitdata, predictdata, p, q, x, n)
    }
  }
  
  ffamily <- f(object = object, fitdata = fitdata, predictdat = predictdata)
  
  
  if(censtype == 'missing'){
    return(list(
      draw = ffamily(func = 'r', n = nrow(predictdata)),
      quantiles = apply(quantprob,
                        MARGIN = 2,
                        FUN = function(q) ffamily(func = 'q', p = q)))
      )
    
  } else if (censtype == 'right') {
    # pindex is the cum. probability up until the censored variable:
    pindex <- ffamily(func = 'p', q = predictdata[[censor]])
    
    # Inverse sampling:
    psample <- runif(n = nrow(predictdata), min = pindex, max = 1)
    draw <- ffamily(func = 'q', p = psample)
    
    # Remap the quantiles on the applicable region in the cdf (psample):
    qindex <- (1-pindex)*quantprob + pindex
    
  }else if (censtype == 'left'){
    # pindex is the cum. probability up until the censored variable.
    # Note that position in psample is reverted to 'right':
    pindex <- ffamily(func =  'p',  q = predictdata[[censor]])
    
    # Inverse sampling:
    psample <- runif(n = nrow(predictdata), min = 0, max = pindex)
    draw <- ffamily(func = 'q', p = psample)
    
    # Remap the quantiles on the applicable region in the cdf (psample):
    qindex <- (pindex-0) * quantprob 
    
  }else if (censtype == 'interval'){ 
    # Applicable inverse sampling region is within the interval:
    pindexupper = ffamily(func = 'p', q = predictdata[[censor]])
    pindexlower = ffamily(func = 'p', q = predictdata[[intervalstart]])
    
    # Inverse sampling:
    psample <- runif(n = nrow(predictdata), min = pindexlower, max = pindexupper)
    draw <- ffamily(func =  'q', p = psample)
    
    # Remap the quantiles on the applicable region in the cdf (psample):
    qindex <- (pindexupper-pindexlower)*quantprob + pindexlower
    
  }else{
    stop('Invalid censtype: Must be one of missing/right/left/interval!')
  }
  
  # Evaluate the rescaled quantiles on the (known) parametrized distribution:
  quantiles <- apply(qindex,
                     MARGIN = 2,
                     FUN = function(q) ffamily(func = 'q', p = q) #  bottleneck?
  )
  
  return(list(
    draw = draw,
    quantiles = quantiles
  ))
}



