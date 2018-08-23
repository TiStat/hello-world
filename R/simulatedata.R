#' @title Simluate missing/censored data
#' @description Data generator for missing/censored data with Normal distribution
#' @param n number of generated observations
#' @param param.formula
#' @param variablenames vector filled with characters specifying all variable
#'  names, that are used in param.formulas. First string is the variable that is
#'  to be censored
#' @param family character. Specifies the gamlss family, from which the data is
#'  drawn. e.g. 'NO' for a dependent variable drawn

#' @param defect formula. specifies which data subset is defected and how.
#'  Formula is specified as follows: ~subset condition | probability | defect.
#'  Subset condition describes which part of the uniformly distributed covariate
#'  [0,1] (whose name must be specified in variablenames in order to be
#'  generated) is selected to be exposed to defection by a bernoulli variable
#'  with the specified probability. The defection funciton is left at the users
#'  convenience. The defect is some function, that the user must specify. For
#'  missing data specifiy NA.
#' @param correlation matrix. If a correlation/covariance matrix is provided, the

#'  covariates specified in variablenames are inversely generated with a
#'  surrogate Multivariate normal distribution to establish correlation between
#'  draws. The draws' cumulative normal probabilities are evaluated in the
#'  Uniform distribution to arrive at correlated uniformly distributed
#'  covariates.
#'  
#' @return List of Dataframes. 'truedata' and 'defected' are dataframes
#'   containing the dependent (generated according to the param.formula list),
#'   the generated covariates, and a censoring/missing 'indicator' The mere
#'   difference between the two Dataframes is, that 'defected' has arteficially generated
#'   censored/missing values according to the 'defect' specification.
#' @examples 
# missing: defect = list(subset = ~ x1 > 0.6, prob = 0.8 , damage = NA)
# right: defect = list(subset = ~ x1 > 0.6, prob = 0.8 , damage = ~1/3*x1)
# rightRandom: defect = list(subset = ~ x1 > 0.6, prob = 0.8 , damage = c(0.01,1))
# left: defect = list(subset = ~ x1 > 0.6, prob = 0.8 , damage = 4/3)
# intervalRandom:defect = list(subset = ~ x1 > 0.6, prob = 0.8 , damage = list(lower=c(0.01,1), upper=c(1.01, 2)))
# 
# subset is valid if it references only the to be defected. this may be if the level of the variable itself 
# changes the probability to be defected. does not alter the independence assumption.
#'@export 

simulateData = function(n,
                        param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)), 
                        defect = list(name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage =1/3),
                        family = 'NO',
                        correlation = NULL) {
  
  # extract all variable names from user specified formula
  varnames = unique(unlist(sapply(param.formula, FUN = all.vars)))
  
  if(!(defect$name %in% varnames)){
    stop('defected variable is not used in param.formula.')
  }
  
  if(!identical(setdiff(names(defect), c('name','subset', 'prob', 'damage')), character(0))){
    stop('defect is not correctly specified')
  }
  
  if(defect$prob>1 | defect$prob<0){
    stop('defect$prob is not a probability.')
  }
  
  if(any(!sapply(param.formula, FUN = function(x) class(x) == 'formula'))){
    stop('at least one param.formula member is not a formula')
  }
  
  blankdata = generateblankdata(varnames, n, correlation)
  
  # evaluate the formulas on data.frame
  param.frame = lapply(param.formula, FUN = function(x) eval(x[[2]], envir = blankdata))
  param.frame$n = n
  
  if(!is.null(param.frame$sigma) && param.frame$sigma<0){
    stop('sigma formula does not ensure positive sigma on all possible covariate values, which are from interval [0,1]')
  } 
  
  # draw TRUE data from conditional family 
  rfam = paste('r', family, sep= '')
  if (any(!names(param.frame) %in% names(formals(rfam)))) {
    stop('the specified parameter formulas do not match the required family\'s parameters')
  }
  y = do.call(rfam, param.frame)
  
  
  defectdata = interpretdamage(blankdata, rule = defect)
  indicator = defectdata$indicator
  return(list(truedata = data.frame(y,blankdata,indicator), 
              defected = cbind(y, defectdata$defected, indicator),
              censtype = defectdata$censtype))
}


#' @description Takes the defect argument and interprets it on the
#' @param rule defect argument of simulate data
#' @return list. Defected dataframe and censtype.
interpretdamage = function(truedata, rule){
  
  # @description function differentiates, if input is a skalar or vektor
  # @return if x is skalar, return x. if x is a vector of length two,
  # drawn n dimensional vector from unifom on interval [x[1], x[2]]
  f = function(x, n){ # differentiate random /fix
    if (length(x) == 1) { # fix
      return(multiply = x)
    } else if (length(x) == 2 ){ # random
      return(multiply = runif(n, x[1], x[2]))
    } else{
      stop('defect$damage length is incorrect')
    }
  }
  
  # subset of observations to be defected with prob.
  indicator = rep(0,nrow(truedata))
  indicator[eval(rule$subset[[2]], envir = truedata)] = 
    rbinom(n = sum(eval(rule$subset[[2]], envir = truedata)), # sum(BOOLEAN)
           size = 1, 
           prob = rule$prob) 
  
  if (all(is.na(rule$damage))) {
    censtype = 'missing'
    truedata[rule$name][indicator == 1,] = NA
    return (list( defected = truedata, 
                  censtype = censtype,
                  indicator = indicator))
  }
  
  n = sum(indicator) # sum(boolean)
  subsetdefect =  truedata[rule$name][indicator == 1,]
  
  if (is.list(rule$damage)) {
    censtype = 'interval'
    truedata$lower = 0
    truedata$lower[indicator == 1] = f(rule$damage[[1]], n) * subsetdefect
    truedata[rule$name][indicator == 1,] = f(rule$damage[[2]], n) * subsetdefect
    
  } else {
    multiply = f(rule$damage, n)
    censtype = ifelse(all(multiply < 1), 'right', 'left')
    truedata[rule$name][indicator == 1,] = multiply * subsetdefect
  }
  return (list( defected = truedata, 
                censtype = censtype,
                indicator = indicator))
}

#' @description This function generates a blank data.frame with uniformly
#'   distributed variables specified in varnames. The variables are drawn i.i.d.
#'   uniform distrib. if correlation = NULL, and are drawn from a correlated
#'   uniform distrib. if a correlation matrix is supplied. The algorithm for the
#'   introduction of correlation is follows to:
#'   https://www.r-bloggers.com/easily-generate-correlated-variables-from-any-distribution-without-copulas/
#'
#' @param n integer. Number of observations
#' @param varnames character vector. specifies variables to be created
#' @param correlation symmetric correlation matix
#' @return data.frame
generateblankdata <- function(varnames, n, correlation= NULL) {
  if(!is.null(correlation)){
    
    # check valid correlation matrix
    if(!is.matrix(correlation)){
      stop('correlation is not a matrix')
    } else if(nrow(correlation)!= ncol(correlation)){
      stop('correlation matrix is not square')
    } else if (all(correlation != t(correlation))){
      stop('correlation matrix is not symmetric')
    }
    
    mu = rep(0,length(varnames))
    rvars = mvrnorm(n = n, mu = mu, Sigma = correlation)
    pvars = pnorm(rvars)
    rawdata = data.frame(qunif(pvars))
    
  } else {
    # generate random data with no correlation.
    rawdata = data.frame(matrix(runif(n*length(varnames)), 
                                ncol = length(varnames)))
  }
  names(rawdata) = varnames
  
  return(rawdata)
}

# (working calls)----------------------------------
# 
rinterval = simulateData(n= 300,
                         param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)),
                         defect = list(name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage =list(c(0.3, 0.99), c(1.2,1.5))),
                         family = 'NO',
                         correlation = NULL)
# 
# rright = simulateData(n= 300,
#                       param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)),
#                       defect = list(name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage =c(0.3, 0.9)),
#                       family = 'NO',
#                       correlation = NULL)
# 
# rleft = simulateData(n= 300,
#                      param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)),
#                      defect = list(name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage =c(1.2, 1.5)),
#                      family = 'NO',
#                      correlation = NULL)
# 
# finterval = simulateData(n= 300,
#                          param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)),
#                          defect = list(name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage =list(1/3,4/3)),
#                          family = 'NO',
#                          correlation = NULL)

# fright = simulateData(n= 300,
#                       param.formula = list(mu = ~exp(x1)+ x2+ x3, sigma = ~sin(x2)),
#                       defect = list(name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage =1/3),
#                       family = 'NO',
#                       correlation = matrix(c(1,0.3,0.2,
#                                              0.3,1, 0.4,
#                                              0.2,0.4,1), nrow = 3))

# fleft = simulateData(n= 300,
#                      param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)),
#                      defect = list(name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage =4/3),
#                      family = 'NO',
#                      correlation = NULL)

# missing = simulateData(n= 300,
#                        param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)), 
#                        defect = list(name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage =NA),
#                        family = 'NO',
#                        correlation = NULL)
# any(is.na(missing$defected$x1))





