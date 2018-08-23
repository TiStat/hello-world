#' @title Imputing censored covariates - GAMLSS

#' @description The MICE Algorithm (Multiple Imputation by Chained Equations) is
#'   a method to impute missing data. This function uses this algorithm for
#'   imputing censored data, using inverse sampling to utilize the additional
#'   information.
#' @param data data.frame containing a dummy censoring indicator, 0 if not
#'   censored/missing, 1 if censored/missing. Note that in case of right (left)
#'   censoring, the censored variable contains the respective minmal (maximal)
#'   duration of the followup which is used for conditional imputation. In case
#'   of interval censored data, two columns are required; specifying the start
#'   and end durations of the interval in question. This implementation assumes,
#'   that the start duration is the observed time, in which no failure occured
#'   before the interval is entered at which the exact point of the observed
#'   state change is unknown. For inverse sampling, the distribution is
#'   conditoned on the start duration and cut at the end duration to ensure the
#'   constrained.
#' @param indicator character. Name of dummy column in data, which indicates the
#'   damaged observation.
#' @param intervalstart character. Name of the column of interval starting
#'   values. by convention, the starting duration in this column is assumed to
#'   be the time passed without failure, before entering the interval, in which
#'   the exact time of failure is unknown.
#' @param censtype character. type of the damaged observation 'missing',
#'   'right', 'left' or 'interval'
#' @param xmu_formula formula for location parameter of xfamily. Dependent
#'   variable specifys the variable which is partially censored/missing and is
#'   to be imputed
#' @param xsigma_formula formula 
#' @param xnu_formula formula 
#' @param xtau_formula formula  
#' @param xfamily gamlss family object.
#' @param ... additional arguments passed in all gamlss fit.
#' 
#' @return 
#' @example 
#' missing = simulateData(n = 100, param.formula = list(mu = ~exp(x1),
#' sigma = ~sin(x2)), variablenames =  c('x1', 'x2'), defect = ~ x1 > 0.6 | 0.8
#' | NA, family = 'NO') 
#' imputex(data = missing$defected, xmu_formula= x1~y+x2,
#' xsigma_formula = ~1, xnu_formula = ~1, xtau_formula = ~1, xfamily =
#' NO(mu.link = 'identity'), indicator = "indicator", censtype= 'missing' )
#' @export
imputex <- function(xmu_formula,
                    xsigma_formula = ~1,
                    xnu_formula = ~1,
                    xtau_formula = ~1,
                    xfamily = NO(mu.link = 'identity'),
                    data,
                    indicator,
                    censtype = c('missing', 'right', 'left', 'interval'),
                    intervalstart = NULL,
                    ...){
  censtype <- match.arg(censtype)
  
  if(!(is.data.frame(data) && !nrow(data) == 0)){
    stop('data must be (non empty) data.frame')
  }
  if(length(xmu_formula)!= 3){
    stop('xmu_formula must be specified as: censoredcovariate ~ set Of covariates')
  }
  
  if(!(is.character(indicator) && indicator %in% names(data))){
    stop('indicator must be a column name in data')
  }
  
  if(!censtype == 'interval' && !is.null(intervalstart) ){
    stop('intervalstart is not required for estimation')
  } else if(censtype == 'interval' && is.null(intervalstart)){
    stop('intervalstart must be specified')
  }
  
  # split data set. for more ambitious projects, this function must be remastered
  # to cope with multiple imputations. In this case, a clear hirarchy of subsetting
  # and imputing should be established. NOTE that this even further increases 
  # the uncertainty contained in the full dataset. 
  W <- function(data, indicator){ 
    df_obs <- data[data[indicator] == 0, ]
    df_cens <- data[data[indicator] == 1, ]
    return(list(obs = df_obs, cens = df_cens))
  }
  
  # split dataset in fully observed & missing/censored data
  Wdat <- W(data, indicator)
  censor <- as.character(xmu_formula[[2]])
  
  # Step 1: fit gamlss with user specified xfamily and formula on observed data
  obsmodel <- gamlss(
    formula = xmu_formula,
    sigma_formula = xsigma_formula,
    nu_formula = xnu_formula,
    tau_formula = xtau_formula,
    family = xfamily,
    data = Wdat$obs,
    ...)
  
  # Step 2: Resampling from fitted model
  # note that these are independent draws from the same distribution. Reframe to m vectors:
  draws <- family_fun(object = obsmodel,
                      func = "r",
                      fitdata = Wdat$obs,
                      predictdata = Wdat$obs,
                      n = nrow(Wdat$obs)*nrow(Wdat$cens))
  
  draws <- data.frame(matrix(draws,
                             nrow = nrow(Wdat$obs),
                             ncol = nrow(Wdat$cens)))
  
  # Basis for respective Bootstap samples (draws[,j], W).
  drawsW <- cbind(draws, Wdat$obs)
  
  # Bootstap samples on simulated observations (rowwise)
  boot <- drawsW[sample(x= 1:nrow(drawsW),
                        size = nrow(Wdat$obs),
                        replace = TRUE), ]
  
  # step 3 estimate param. on each respective set {x*boot(j), W_obs} for all j
  bootmodel <- list()
  imputemat <- data.frame(X1 = vector(length= nrow(Wdat$cens)))
  imputeq <- list()
  quantil <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  for (i in 1:ncol(draws)){
    # iterate only over the names of booted vectors.
    
    # manipulate the formula for each estimation to get the respective
    # booted x column regressed on W
    bootformula <- as.formula(paste(names(boot)[i], '~',
                                    as.character(as.vector(xmu_formula)[3]),
                                    sep = ''))
    
    bootmodel[[i]] <- gamlss(formula = bootformula,
                             sigma_formula = xsigma_formula,
                             nu_formula = xnu_formula,
                             xtau_formula = xtau_formula,
                             family = xfamily,
                             data = boot,
                             ...)
    
    # Simulate data from the corresponding fitted distribution.
    imputecandidate <- names(boot)[i]
    
    impute <- samplecensored(object = bootmodel[[i]],
                            censtype,
                            fitdata = boot,
                            predictdata = Wdat$cens,
                            censor, # censor becomes upper bound if !is.null(intervallstart)
                            intervalstart = intervalstart,
                            quantil)
    
    
    imputemat[[imputecandidate]] <- impute$draw
    imputeq[[i]] <- impute$quantiles
  }
  # imputed vector
  imputedx <- apply(imputemat, MARGIN = 1,median)
  
  # (output augmentation)-------------------------------------------------------
  # save censored values before they get overwritten
  censoredx <- Wdat$cens[[censor]]
  
  # distances from imputed to censored values
  distances <- switch((censtype == "left" | censtype == "right") + 1,
                      NULL,
                      abs(censoredx - imputedx))
  
  # complete data with imputations
  Wdat$cens[censor] <- imputedx
  fulldata <- rbind(Wdat$obs,Wdat$cens)
  
  # variability of imputed observation among all drawn from booted
  imputevariance <- apply(imputemat, MARGIN = 1, FUN = var)
  
  # average imputed quantiles
  A <- array(unlist(imputeq), 
            dim = c(nrow(imputeq[[1]]),
                    ncol(imputeq[[1]]), 
                    length(imputeq)))
  impquantiles <- as.data.frame(apply(A, c(1,2), mean))
  colnames(impquantiles) <- c('q5','q25','q50', 'q75', 'q95' )
  
  mcall <- match.call()
  
  result <- list(imputations = imputemat,
                 imputedx = imputedx,
                 fulldata = fulldata,
                 mcall = mcall,
                 nimputations = nrow(Wdat$cens),
                 nobservations = nrow(data),
                 censtype = censtype,
                 imputevariance = imputevariance,
                 impquantiles = impquantiles,
                 distances = distances
                 )
  
  #  Create a class for this kind of result
  class(result) <- "imputed"
  
  return(result)
}


