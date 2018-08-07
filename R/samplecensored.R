#' Title Evaluate gamlss family functions from gamlss object
#' @description Function that, given a "gamlss" object, evaluates the
#'   distribution under predicted parameters of the provided dataframe
#'   "predictdata".
#'   CAUTION: Exactly ONE of the arguments x,q,p,n MUST be specfied!
#' @param object gamlss fit object
#' @param func character. "d", "p", "q", "r" for either density, distr. function,
#'   quantile or random generation
#' @param fitdata
#' @param predictdata dataframe. Containing the observations for which the
#'   parameters are predicted.
#' @param x,q scalar numeric. Quantile value. If density or probability function
#'   used respectively
#' @param p scalar numeric. Probability value. If quantile function used
#' @param n scalar numeric. Number of observations. If random generator function used
#' @param ... argumenst to be passed to the called distributional function
#'
#' @return Depending on the choice of func, the respective vector of
#'   (d)density-, (p)cdf.- , (q)quantile- or (r)random values is returned
family_fun <- function(object, func, fitdata, predictdata ,p = NULL, q = NULL, x = NULL, n = NULL, ...) {

  # find the correct family function to evaluate
  fam_name = object$family[1]
  f_fun = paste(func, fam_name, sep= '')

  prd <- predictAll(
    object = object,
    type = "response",
    newdata = predictdata,
    data = fitdata
  )

  param = list(
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
  param = param[!sapply(param, is.null)] # kick out NULL parameters

  if (any(!names(param) %in% names(formals(f_fun)))) {
    stop("One of x,q,p,n,... arguments doesn't match with the distributional function
         (e.g. dNO, pNO, qNO, rNO). See the family's documentation for admissable arguments.")
  }
  return(do.call(f_fun, param))
  }

# Beispiel für die error message
# family_fun(gamlssF, func = 'q',predictdata = predict.df, x = 0)
# family_fun(gamlssF, func = 'q',predictdata = predict.df, p = 0.5)


#' @description Inverse sampling of censored variables, to impute only valid observations
#' conditional on the respective fit
#' @param object gamlss. Fitted model whose parameters are predicted for predictdata.
#' @param censtype chr. specifies the type of censoring, for which is imputed
#' @param predictdata dataframe. predictdata of the censored observations, for which
#'   imputations are drawn
#' @param fitdata dataframe. The orignal Dataset upon which gamlss was fitted
#' @param censor character. Name of the to be predicted (damaged) column in predictdata.
#' is only required if censtype is NOT 'missing'.
#' @return
samplecensored = function(object, censtype, predictdata, fitdata, censor){  # predictdata is W$obs/W$mis, an denen die fitted values predicted werden, um anschließend p auszuwerten
  if(censtype == 'missing'){
    return(family_fun(object, func = 'r', predictdata = predictdata, n = nrow(predictdata), fitdata = fitdata))

  }else if(censtype == 'right'){
    pindex = family_fun(object, func='p', predictdata = predictdata, q= predictdata[[censor]],  fitdata = fitdata)
    psample = runif(n = nrow(predictdata), min = pindex, max = 1)
    return(family_fun(object, func =  'q', predictdata = predictdata, p = psample,  fitdata = fitdata))

  }else if (censtype == 'left'){
    pindex = family_fun(object, func = 'p', predictdata = predictdata , q=predictdata[[censor]],  fitdata = fitdata)
    psample = runif(n = nrow(predictdata), min = 0, max = pindex)
    return(family_fun(object, func = 'q',predictdata = predictdata, p = psample,  fitdata = fitdata))

  # }else if (censtype == 'interval'){ input format noch unklar

  }else{
    stop('invalid method')
  }
}



