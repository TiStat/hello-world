#' @title Simluate missing/censored data
#' @description Data generator for missing/censored data with Normal distribution
#' @param n number of generated observations
#' @param param.formula
#' @param variablenames vector filled with characters specifying all variable
#'  names, that are used in param.formulas. First string is the variable that is
#'  to be censored
#' @param family character. Specifies the gamlss family, from which the data is
#'  drawn. e.g. 'NO' for a dependent variable drawn
#' @param name character. specifies variable name to be defected.
#' @param subset formula. states a condition ( e.g. ~x1 > 0.6) which specifies
#'   the fraction of observations, that are to be defected. Note, that if
#'   'subset' does not exclusivly use the 'name[d]' variable, this implies that
#'   the independence assumption of MICE is not met (on purpose). e.g. of unmet
#'   condition: (x2 < 0.3 & x3<0.2)
#' @param prob numeric value. Specifies the binomial probability for each
#'   observation in 'subset' to be defected.
#' @param damage By users defintion, it specifies what type and how the data is
#'   to be defected.
#'   'damage' = NA generates missing data. 
#'   A value between [0, 1] implies right censoring (e.g. 'damage' = 1/3), 
#'   [1,...] left censoring. The value is used to multiply the truevalue of 
#'   'name' in order to defect the data. 
#'   The generalization for fixed interval factors is 'damage' = list(1/3, 4/3),
#'   where the values specifiy the factor for the lower and the upper bound respectively.
#'   More realistic examples can be generated with vector valued 'damage':
#'   If 'damage' = c(0.1, 1) is a vector of length 2, it specifies the min 
#'   and max value of a uniform distribution, from which a factor is randomly 
#'   drawn for each observation with which the true data is multiplied.
#'   The generalization for random interval factors is 'damage' = list(c(0.2, 1), c(1,3)),
#'   where the first vector specifies the unif interval for factors affecting 
#'   the lower bound and the second affecting the upper bound.
#'   NOTE: if a list is provided, both members must either vectors or sigle values.
#' @param correlation matrix. If a correlation/covariance matrix is provided,
#'   the drawn variables are uniformly drawn, but correlated according to this matrix. 
#'  
#' @return List of Dataframes. 'truedata' and 'defected' are dataframes
#'   containing the dependent (generated according to the param.formula list),
#'   the generated covariates, and a censoring/missing 'indicator' The mere
#'   difference between the two Dataframes is, that 'defected' has arteficially generated
#'   censored/missing values according to the 'defect' specification.
#' @examples 
# missing: damage = NA
# right: damage = ~1/3*x1
# rightRandom:  damage = c(0.01,1)
# left:  damage = 4/3
# intervalfix: damage = list(1/3, 4/3)
# intervalRandom: damage = list(c(0.01,1), c(1.01, 2))
#'@export 

simulateData = function(n,
                        param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)), 
                        name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage =1/3,
                        family = 'NO',
                        correlation = NULL) {
  
  # extract all variable names from user specified formula
  varnames = unique(unlist(sapply(param.formula, FUN = all.vars)))
  
  if(!(name %in% varnames)){
    stop('defected variable is not used in param.formula.')
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
  
  
  defectdata = interpretdamage(blankdata, name, subset, prob, damage)
  indicator = defectdata$indicator
  return(list(truedata = data.frame(y,blankdata,indicator), 
              defected = cbind(y, defectdata$defected, indicator),
              censtype = defectdata$censtype))
}

#' @title Defect existing data.
#' @description Takes the defect argument and interprets it on the truedataset.
#' Note, that interpretdamage is generic and can defect any existing dataset according to rule
#' @param name character. specifies variable name to be defected.
#' @param subset formula. states a condition ( e.g. ~x1 > 0.6) which specifies
#'   the fraction of observations, that are to be defected. Note, that if
#'   'subset' does not exclusivly use the 'name[d]' variable, this implies that
#'   the independence assumption of MICE is not met (on purpose).
#' @param prob numeric value. Specifies the binomial probability for each
#'   observation in 'subset' to be defected.
#' @param damage By users defintion, it specifies what type and how the data is
#'   to be defected.
#'   'damage' = NA generates missing data. 
#'   A value between [0, 1] implies right censoring (e.g. 'damage' = 1/3), 
#'   [1,...] left censoring. The value is used to multiply the truevalue of 
#'   'name' in order to defect the data. 
#'   The generalization for fixed interval factors is 'damage' = list(1/3, 4/3),
#'   where the values specifiy the factor for the lower and the upper bound respectively.
#'   More realistic examples can be generated with vector valued 'damage':
#'   If 'damage' = c(0.1, 1) is a vector of length 2, it specifies the min 
#'   and max value of a uniform distribution, from which a factor is randomly 
#'   drawn for each observation with which the true data is multiplied.
#'   The generalization for random interval factors is 'damage' = list(c(0.2, 1), c(1,3)),
#'   where the first vector specifies the unif interval for factors affecting 
#'   the lower bound and the second affecting the upper bound.
#'   NOTE: if a list is provided, both members must either vectors or sigle values.
#' @return list. Defected dataframe, an indicator vector specifying which
#'   observation was defected and atomatically description of censtype, which
#'   the user implicitly defined by damage.
#' @export 
interpretdamage = function(truedata, name, subset, prob, damage){
  
  # catch corner cases of user failure owed to complex input.
  if(class(subset) != 'formula' ){
    stop('formula must be a formula')
  }
  if(prob>1 | prob<0){
    stop('prob is not a probability.')
  }
  
  if(is.list(damage)){ # interval case
    if(length(unique(sapply(damage, FUN = length)))!= 1){
      stop('both damage members must be of same length: either vector or single values')
    } else if (any(sapply(damage, FUN = length) > 2)) {
      stop('the vectors in damage list must both be either of length 1 or 2')
    } 
    
    # die uniform bounds  <1 >1 im interval case
    if(!all(all(damage[[1]]<1), all(damage[[1]]>0))){
      stop('lower bound interval factor is not from interval (0,1) excluding {0,1} ')
    }
    if(!all(damage[[2]]>1)){
      stop('upper bound interval factor is not from interval (1,INF) excluding {1} ')
    }
  }

  #' @description function differentiates, if input is a skalar or vektor
  #' @return if x is skalar, return x. if x is a vector of length two, drawn n
  #'   dimensional vector from unifom on interval [x[1], x[2]]
  f = function(x, n){ # differentiate random /fix
    if (length(x) == 1) { # fix
      return(multiply = x)
    } else if (length(x) == 2 ){ # random
      return(multiply = runif(n, x[1], x[2]))
    }
  }
  
  # subset of observations to be defected with prob.
  indicator = rep(0,nrow(truedata))
  indicator[eval(subset[[2]], envir = truedata)] = 
    rbinom(n = sum(eval(subset[[2]], envir = truedata)), # sum(BOOLEAN)
           size = 1, 
           prob = prob) 
  
  if (all(is.na(damage))) {
    censtype = 'missing'
    truedata[name][indicator == 1,] = NA
    return (list( defected = truedata, 
                  censtype = censtype,
                  indicator = indicator))
  }
  
  n = sum(indicator) # sum(boolean)
  subsetdefect =  truedata[name][indicator == 1,]
  
  if (is.list(damage)) {
    censtype = 'interval'
    truedata$lower = 0
    truedata$lower[indicator == 1] = f(damage[[1]], n) * subsetdefect
    truedata[name][indicator == 1,] = f(damage[[2]], n) * subsetdefect
    
  } else {
    multiply = f(damage, n)
    censtype = ifelse(all(multiply < 1), 'right', 'left')
    truedata[name][indicator == 1,] = multiply * subsetdefect
  }
  return (list( defected = truedata, 
                censtype = censtype,
                indicator = indicator))
}


#' @title Generating a blank data frame
#' @description This function generates a blank data.frame with uniformly
#'   distributed variables specified in varnames. The variables are drawn i.i.d.
#'   uniform distrib. if correlation = NULL, and are drawn from a correlated
#'   uniform distrib. if a correlation matrix is supplied. The algorithm for the
#'   introduction of correlation is follows to:
#'   https://www.r-bloggers.com/easily-generate-correlated-variables-from-any-distribution-without-copulas/
#'   The idea is: covariates specified in variablenames are inversely generated
#'   with a surrogate Multivariate normal distribution to establish correlation
#'   between mvnormal draws. The draws' cumulative normal probabilities are
#'   evaluated in the Uniform distribution to arrive at correlated uniformly
#'   distributed covariates.
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
    rvars = MASS::mvrnorm(n = n, mu = mu, Sigma = correlation)
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





