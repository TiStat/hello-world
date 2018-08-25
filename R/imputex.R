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
#' @param m Number of imputations. Default is m = 5
#' @return Returns internal results of the algorithm.
#' @examples 
#' missing = simulateData(n = 100, param.formula = list(mu = ~exp(x1) + x3,
#' sigma = ~sin(x2)), name = 'x1', subset = ~ x1 > 0.6, prob = 0.8,
#' damage = NA, family = 'NO', correlation = NULL) 
#' imputex(data = missing$defected, xmu_formula= x1~y+x3,
#' xsigma_formula = ~x2, xnu_formula = ~1, xtau_formula = ~1, xfamily =
#' NO(mu.link = 'identity'), indicator = "indicator", censtype= 'missing')
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
                    m = 5,
                    ...){
  censtype <- match.arg(censtype)
  
  if(!(is.data.frame(data) & !nrow(data) == 0)){
    stop('data must be (non empty) data.frame')
  }
  if(length(xmu_formula)!= 3){
    stop('xmu_formula must be specified as: censoredcovariate ~ set Of covariates')
  }
  
  if(!(is.character(indicator) & indicator %in% names(data))){
    stop('indicator must be a column name in data')
  }
  
  if(!censtype == 'interval' & !is.null(intervalstart) ){
    stop('intervalstart is not required for estimation')
  } else if(censtype == 'interval' & is.null(intervalstart)){
    stop('intervalstart must be specified')
  }
  

  # split dataset in fully observed & missing/censored data
  Wdat <- list(obs = data[data[indicator] == 0, ],
               cens = data[data[indicator] == 1, ])
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
  
  # Step 2: Resampling from fitted model.
  # Note that these are m independent draws from vectors of length nrow(Wdat$obs), which are stacked!
  # Reframe to m vectors:
  draws <- family_fun(object = obsmodel,
                      func = "r",
                      fitdata = Wdat$obs,
                      predictdata = Wdat$obs,
                      n = nrow(Wdat$obs) * m)
  
  # Be careful when unstacking drawn vectors!
  draws <- data.frame(matrix(draws,
                             nrow = nrow(Wdat$obs),
                             ncol = m))
  
  # Basis for respective Bootstap samples (draws[,j], W).
  drawsW <- cbind(draws, Wdat$obs)
  
  # Bootstap samples on simulated observations (rowwise)
  boot <- drawsW[sample(x= 1:nrow(drawsW),
                        size = nrow(Wdat$obs),
                        replace = TRUE), ]
  
  # step 3 estimate param. on each respective set {x*boot(j), W_obs} for all j
  bootmodel <- list()
  proposals <- data.frame(X1 = vector(length= nrow(Wdat$cens)))
  imputeq <- list()
  quantil <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  for (i in 1:m){
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
    
    
    proposals[[imputecandidate]] <- impute$draw
    imputeq[[i]] <- impute$quantiles
  }
  # imputed vector
  imputedx <- apply(proposals, MARGIN = 1,median)
  
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
  imputevariance <- apply(proposals, MARGIN = 1, FUN = var)
  
  # average imputed quantiles
  A <- array(unlist(imputeq), 
            dim = c(nrow(imputeq[[1]]),
                    ncol(imputeq[[1]]), 
                    length(imputeq)))
  imputequantiles <- as.data.frame(apply(A, c(1,2), mean))
  colnames(imputequantiles) <- c('q5','q25','q50', 'q75', 'q95' )
  
  mcall <- match.call()
  
  result <- list(proposals = proposals,
                 imputedx = imputedx,
                 fulldata = fulldata,
                 mcall = mcall,
                 nobservations = nrow(data),
                 nreplacements = nrow(Wdat$cens),
                 censtype = censtype,
                 imputevariance = imputevariance,
                 imputequantiles = imputequantiles,
                 distances = distances,
                 m = m,
                 Wobs = Wdat$obs
                 )
  
  #  Create a class for this kind of result
  class(result) <- "imputed"
  
  return(result)
}


