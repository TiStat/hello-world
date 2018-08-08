#'@description Data generator for missing/censored data with Normal distribution
#'@param n number of generated observations
#'@param ymu.formula
#'@param ysigma.formula
#'
#'@param family character. Specifies the gamlss family, from which the data is drawn. e.g. 'NO' for
#'a dependent variable drawn
#'@return Dataframe, containing covariates, censor indicator and the logit value
#'pi which is modeled as surrogate for binary decision on Breastcancer
simulateData = function(method,
                        n,
                        ymu.formula = NULL,
                        ysigma.formula = NULL,
                        ynu.formula = NULL,
                        ytau.formula = NULL,
                        # formula to defect the data
                        family) {

  # True process
  x = runif(n)
  x1 = x
  x2 = runif(n)

  # draw y from family with conditional parameters
  rfam = paste('r', family, sep = '')
  param = list(
    mu = eval(ymu.formula[[2]]),
    sigma = eval(ysigma.formula[[2]]),
    nu = eval(ynu.formula[[2]]),
    tau = eval(ytau.formula[[2]]),
    n
  )
  param = param[!sapply(param, is.null)]
  y = do.call(rfam, param)

  indicator = rbinom(n,1, prob= 0.05)

  # alter data according to defect method
  if (method == 'missing'){
    x1[indicator == 1] = NA
  }else if (method == 'right'){
    x1[indicator == 1] =  1/3 * x1[indicator == 1]
  }else if (method == 'left'){

  }else if (method == 'interval'){

  }else{
    stop('invalid method is specified')
  }
  return(list(defected = data.frame(y,x1,x2, indicator), truedata = data.frame(y,x,x2, indicator)))

}

# simulated package run
library(ggplot2)
run <- function(censtype,
                indicator,
                n,
                ymu.formula,
                ysigma.formula,
                family = 'NO') {

  data = simulateData(method = censtype ,
                      n,
                      ymu.formula,
                      ysigma.formula,
                      family = 'NO')

  d2 <- imputex(data = data$defected,
                xmu_formula= x1~y+x2,
                xsigma_formula = ~1,
                xnu_formula = ~1,
                xtau_formula = ~1,
                xfamily = NO(mu.link = 'identity'),
                indicator = "indicator",
                censtype )

  # remaster: dataframe usage is not parsimonious
  # carefull if user specifies more than one dimension.
  p = ggplot() +
    geom_point(data = d2$fulldata, aes(y = y , x = x1, color = indicator)) +
    geom_point(data = data$defected[data$defected[indicator] == 1, ],
               aes(y = y, x = x1),
               color = 'red') +
    geom_point(data = data$truedata[data$truedata[indicator] == 1, ],
               aes(y = y, x = x), color = 'green')

  print(p)

  # fity = fitfull(ymu_formula = y~x1, # hard coded!
  #                ysigma_formula = y~x2,
  #                yfamily= NO(link= 'identity'),
  #                xmu_formula = x1~y,
  #                xfamily = NO(link = 'identity'),
  #                data = d2$fulldata,
  #                indicator = 'inticator',
  #                censtype = 'right')
  return(list(data = data, imputex =  d2)) #, fit = fity ))
}

r = run(
  censtype = 'right',
  indicator = 'indicator',
  n = 1000,
  ymu.formula = ~ 2 * x1 ,
  ysigma.formula =  ~ 0.5 * x2
)
