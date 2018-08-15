#'@description Data generator for missing/censored data with Normal distribution
#'@param n number of generated observations
#'@param param.formula list
#'@param variablenames vector filled with characters specifying all variable names, that are used in 
#'param.formulas. First string is the variable that is to be censored
#'
#'@param family character. Specifies the gamlss family, from which the data is drawn. e.g. 'NO' for
#'a dependent variable drawn
#'@param defect formula. specifies which data subset and how it is defected.
#'is specified as follows: ~subset condition | probability | defect. Subset condition describes which part of the uniformly distributed covariate [0,1] (whose name must be specified in variablenames)
#'is selected to be exposed to defection by a bernoulli variable with the specified probability. The defection funciton is left at the users convenience
#'The defect is some function, that the user must specify. For missing specifiy NA.
#' #'@return Dataframe, containing covariates, censor indicator
#' #'@examples defect = ~ x1 > 0.6 | 0.8 | 1/3*x1
#' missing = simulateData(n = 100, param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)),
#' #'                        variablenames =  c('x1', 'x2'), defect = ~ x1 > 0.6 | 0.8 | NA, family = 'NO')
#' right = simulateData(n = 100, param.formula = list(mu = ~exp(x1)+ x2, sigma = ~sin(x2)),
#'                      variablenames =  c('x1', 'x2'), defect = ~ x1 > 0.6 | 0.8 | 1/3*x1, family = 'NO')
#'@export 
simulateData = function(n,
                        param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)), # ensure sigma positive!
                        variablenames = c('x1', 'x2'),
                        defect = ~ x1 > 0.6 | 0.8 | 1/3*x1,
                        family = 'NO') {
  
  if(any(!sapply(param.formula, FUN = function(x) class(x) == 'formula'))){
    stop('at least one param.formula members is not a formula')
  }
  if(class(defect) != 'formula'){
    stop('defect is not a formula')
  }
  
  # generate some random data
  rawdata = data.frame(matrix(runif(n*length(variablenames)), ncol = length(variablenames)))
  names(rawdata) = variablenames
  
  # evaluate the formulas on data.frame
  param.frame = list()
  for (i in 1:length(param.formula)){
    param.frame[[i]] = eval(param.formula[[i]][[2]], envir = rawdata)
  }
  names(param.frame) = names(param.formula)
  param.frame$n = n

  # draw from conditional family 
  rfam = paste('r', family, sep= '')
  if (any(!names(param.frame) %in% names(formals(rfam)))) {
    stop('the specified parameter formulas do not match the required family\'s parameters')
  }
  y = do.call(rfam, param.frame)
  
  # data = data.frame(y, param.frame, indicator = 0)
  prob = defect[[2]][[2]][[3]] # prob condition
  indicator = rep(0,nrow(rawdata))
  
  # figure out which observation to be defected
  indicator[eval(defect[[2]][[2]][[2]], envir = rawdata)] = 
    rbinom(n = sum(eval(defect[[2]][[2]][[2]], envir = rawdata)),1, prob)
  
  newdata = rawdata
  newdata$x1[indicator == 1] = eval(defect[[2]][[3]], envir = rawdata)[indicator == 1]
   
  return(list(defected = data.frame(y,rawdata, indicator), 
              truedata = data.frame(y,newdata, indicator)))

}
missing = simulateData(n = 100, param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)),
                       variablenames =  c('x1', 'x2'), defect = ~ x1 > 0.6 | 0.8 | NA, family = 'NO')
right = simulateData(n = 100, param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)), 
                     variablenames =  c('x1', 'x2'), defect = ~ x1 > 0.6 | 0.8 | 1/3*x1, family = 'NO')

