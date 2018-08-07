#'@description Data generator for censored data with Normal distribution
#'@param n number of generated observations
#'@param mu.formula
#'@param sigma.formula
#'@param family
#'@return Dataframe, containing covariates, censor indicator and the logit value
#'pi which is modeled as surrogate for binary decision on Breastcancer
generateData = function(method ,n, mu.formula = NULL, sigma.formula = NULL, family = NULL){
  # True process
  x = runif(n)
  x1 = x
  x2 = runif(n)
  y = rnorm(n,mean = eval(mu.formula[[2]]), sd = eval(sigma.formula[[2]]))
  indicator = rbinom(n,1, prob= 0.05)


  if (method == 'missing'){
    x1[indicator == 1] = NA
  }else if (method == 'right'){
    x1[indicator == 1] =  1/3 * x1[indicator == 1]
  }else if (method == 'left'){

  }else if (method == 'interval'){

  }else{
    stop('invalid method is specified')
  }
  return(list(defecteddata = data.frame(y,x1,x2, indicator), truedata = data.frame(y,x,x2, indicator)))

}

# simulated package run
library(ggplot2)
run <- function(censtype, indicator, ...) {
  data = generateData(method = censtype , 1000, ~2*x1 ,~0.5*x2, ...)

  d2 <- imputex(data = data$defecteddata,
                xmu_formula= x1~y+x2,
                xsigma_formula = ~1,
                xnu_formula = ~1,
                xtau_formula = ~1,
                xfamily = NO(mu.link = 'identity'),
                indicator = "indicator",
                censtype )

  ggplot()+
    geom_point(data = d2$fulldata, aes(y = y , x = x1, color = indicator))+
    geom_point(data = data$defecteddata[data$defecteddata$indicator == 1,],
               aes(y = y, x = x1), color = 'red')+
    geom_point(data = data$truedata[data$truedata$indicator == 1,],
               aes(y = y, x = x), color = 'green')
}

run(censtype = 'right', indicator = 'indicator')
