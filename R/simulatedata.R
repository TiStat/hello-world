#'@description Data generator for censored data with Normal distribution
#'@param n number of generated observations
#'@param mu.formula
#'@param sigma.formula
#'@param family
#'@return Dataframe, containing covariates, censor indicator and the logit value
#'pi which is modeled as surrogate for binary decision on Breastcancer
generateData = function(method ,n, mu.formula = NULL, sigma.formula = NULL, family = NULL){
  # True process
  x1 = runif(n)
  x2 = runif(n)
  y = rnorm(n,mean = eval(mu.formula[[2]]), sd = eval(sigma.formula[[2]]))
  indicator = rbinom(n,1, prob= 0.05)
  if (method == 'cens'){
    x1[indicator == 1] =  1/3 * x1[indicator == 1]
  } else if (method == 'missing'){
    x1[indicator == 1] = NA
  }
  return(data.frame(y,x1,x2, indicator))

}
data = generateData(method = 'cens', 1000, ~2*x1 ,~0.5*x2)

library(ggplot2)
qplot(data = data, y = y , x = x1, color = indicator)


