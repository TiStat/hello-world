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
  
  if(!identical(setdiff(names(defect), c('name','subset', 'prob', 'damage')), character(0))){
    stop('defect is not correctly specified')
  }
  
  if(defect$prob>1 | defect$prob<0){
    stop('defect$prob is not a probability.')
  }
  
  if(any(!sapply(param.formula, FUN = function(x) class(x) == 'formula'))){
    stop('at least one param.formula member is not a formula')
  }
  
  if (is.list(defect$damage) && !all(sapply(defect$damage, FUN = is.numeric))){
    stop('elements of defect$damage must be of same type')
  }
  # extract all variable names from user specified formula
  varnames = sapply(param.formula, FUN = all.vars)
  
  if(!(defect$name %in% unlist(varnames))){
    stop('defected variable is not used in param.formula.')
  }
  
  if(!is.null(correlation) && 
     nrow(correlation)!= unique(unlist(varnames))&& 
     nrow(correlation)!= ncol(correlation)){
    stop('The specified correlation matrix is either not square, its nrow is larger than the number of specified variables in param.formula')
  }
  
  # draw correlated data 
  # according to: https://www.r-bloggers.com/easily-generate-correlated-variables-from-any-distribution-without-copulas/
  if(!is.null(correlation)){
    mu = rep(0,length(unique(unlist(varnames))))
    rvars = mvrnorm(n = n, mu = mu, Sigma = correlation)
    pvars = pnorm(rvars)
    rawdata = data.frame(qunif(pvars))
    
  } else {
    # generate random data with no correlation.
    rawdata = data.frame(matrix(runif(n*length(varnames)), 
                                ncol = length(varnames)))
  }
  names(rawdata) = unique(unlist(varnames))
  
  # evaluate the formulas on data.frame
  param.frame = list()
  for (i in 1:length(param.formula)){
    param.frame[[i]] = eval(param.formula[[i]][[2]], envir = rawdata)
  }
  names(param.frame) = names(param.formula)
  
  if(!is.null(param.frame$sigma) && param.frame$sigma<0){
    stop('sigma formula does not ensure positive sigma on all possible covariate values, which are uniformly distributed on the interval [0,1]')
  } 
  
  param.frame$n = n
  
  # draw TRUE data from conditional family 
  rfam = paste('r', family, sep= '')
  if (any(!names(param.frame) %in% names(formals(rfam)))) {
    stop('the specified parameter formulas do not match the required family\'s parameters')
  }
  y = do.call(rfam, param.frame)
  
  # subset of observations to be defected with prob.
  indicator = rep(0,nrow(rawdata))
  indicator[eval(defect$subset[[2]], envir = rawdata)] = 
    rbinom(n = sum(eval(defect$subset[[2]], envir = rawdata)),
           size = 1, 
           prob = defect$prob) # note sum(BOOLEAN)
  
  newdata = rawdata
  if(any(is.na(defect$damage))){ # MISSING
    newdata[defect$name][indicator== 1,] = NA
    
  } else if (!is.list(defect$damage) && length(defect$damage)==1){ # LEFT/RIGHT fixed factor
    newdata[defect$name][indicator== 1,] = 
      newdata[defect$name][indicator == 1,]*defect$damage
    
  } else if (!is.list(defect$damage) && length(defect$damage) == 2){ # LEFT/RIGHT random factors
    newdata[defect$name][indicator == 1,] = 
      newdata[defect$name][indicator == 1,] * 
      runif(n = sum(indicator), min = defect$damage[1], max = defect$damage[2])
    
  } else if (is.list(defect$damage)){ # INTERVAL
    newdata$lower = rep(0, n)
    if (all(unlist(lapply(defect$damage, length))==1)){ # fixed factor
      # lower : as START time
      newdata$lower[indicator == 1] = 
        newdata[defect$name][indicator == 1,]*defect$damage[[1]]
      
      # upper : defect$name as END time
      newdata[defect$name][indicator== 1,] = 
        newdata[defect$name][indicator == 1,]*defect$damage[[2]]
      
    }else if (all(unlist(lapply(defect$damage, length))==2)){ # random factor
      # lower
      newdata$lower[indicator == 1] = 
        newdata[defect$name][indicator == 1,] * 
        runif(n = sum(indicator), min = defect$damage[[1]][1], max = defect$damage[[1]][2])
      
      # upper
      newdata[defect$name][indicator == 1,] = 
        newdata[defect$name][indicator == 1,] * 
        runif(n = sum(indicator), min = defect$damage[[2]][1], max = defect$damage[[2]][2])
    }
  }
  return(list(truedata = data.frame(y,rawdata, indicator), 
              defected = data.frame(y,newdata, indicator)))
}


# (working calls)----------------------------------

rinterval = simulateData(n= 300,
                         param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)), 
                         defect = list(name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage =list(c(0.3, 0.99), c(1.2,1.5))),
                         family = 'NO',
                         correlation = NULL)

rright = simulateData(n= 300,
                      param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)), 
                      defect = list(name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage =c(0.3, 0.9)),
                      family = 'NO',
                      correlation = NULL)

rleft = simulateData(n= 300,
                     param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)), 
                     defect = list(name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage =c(1.2, 1.5)),
                     family = 'NO',
                     correlation = NULL)

finterval = simulateData(n= 300,
                         param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)), 
                         defect = list(name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage =list(1/3,4/3)),
                         family = 'NO',
                         correlation = NULL)

fright = simulateData(n= 300,
                      param.formula = list(mu = ~exp(x1)+ x2+ x3, sigma = ~sin(x2)), 
                      defect = list(name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage =1/3),
                      family = 'NO',
                      correlation = matrix(c(1,0.3,0.2,
                                             0.3,1, 0.4,
                                             0.2,0.4,1), nrow = 3))

fleft = simulateData(n= 300,
                     param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)), 
                     defect = list(name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage =4/3),
                     family = 'NO',
                     correlation = NULL)

missing = simulateData(n= 300,
                       param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)), 
                       defect = list(name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage =NA),
                       family = 'NO',
                       correlation = NULL)
any(is.na(missing$defected$x1))



