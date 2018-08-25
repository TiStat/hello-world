#' @title Printing an object of class "imputed"
#' @param x Object of class "imputed"
#'
#' @return Retruns a print in the console
#' @export
print.imputed <- function(x, ...){
  
  v <- x$nreplacements
  cens_type <- x$censtype
  n <- x$nobservations
  mcall <- x$mcall
  
  cat("\n", 
      cat("Call:"),
      paste(deparse(mcall), sep = "\n", collapse = "\n"),
      sprintf("\n\n Number of observations: %s\n", n),
      paste("Replaced", v, cens_type, "censored values") )
}

#' @title Summarizing an object of class "imputed"
#' @param object Object of class "imputed"
#'
#' @return Retruns a summary in the console
#' @export
summary.imputed <- function(object, ...) {
  
  v <- object$nreplacements
  n <- object$nobservations
  sm <- summary(object$imputedx)
  mcall <- object$mcall
  
  cens_type <- object$censtype
  r <- ifelse(!is.null(object$distances), round(mean(object$distances), 3), "undefined")
  
  cat("\n", 
      cat("Call:  \n", paste(deparse(mcall), sep = "\n", collapse = "\n")),
      paste("\n Number of observations:", n),
      paste("\n Type of censoring:", cens_type),
      paste("\n Number of replacements:", v),
      "\n\n Imputed values:", 
      "\n",
      paste(names(sm), collapse = "     "),
      "\n",
      paste(round(sm, 3), collapse = "      "),
      
      paste("\n\n Average distance of imputations to censorings: \n", r)
  )
}


#' @title Plot proposed imputations and their conditional censored
#'   distributions' average quantiles
#' @description Given an imputex object, this funciton plots (1) the the
#'   approximate averaged quantiles of the valid part of the censored conditonal
#'   distributions, from which the proposed vectors were 'drawn', before the
#'   proposals were aggregated via median to become the imputed vector. (2) the
#'   actual imputations for each observation. the red dots are the median
#'   proposals, which are infact the final imputed vector. (3) Density plots on
#'   each Proposalvector. Consider, that each proposalvector is drawn from a
#'   bootmodel-distribution whose observations have their own parameters. So
#'   each proposalvector is a realization of a mixed distribution. Suppose, that
#'   the booted models are relatively simmilar, the density estimates on all
#'   proposalvectors (which are displayed) should be relatively simmilar with
#'   little variation.
#' @param x imputex object.
#' @param boxes boolean. indicating, whether (2) should be displayed as a
#'   boxplot. Note, that the median values are the imputation values.
#' @examples 
#' rinterval = simulateData(n= 300,
#'                          param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)),
#'                          name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 ,
#'                          damage =list(c(0.8, 0.99), c(1.2,1.5)),
#'                          family = 'NO',
#'                          correlation = NULL)
#' d  <- imputex(data = rinterval$defected,
#'               xmu_formula= x1~y,
#'               indicator = "indicator",
#'               censtype = 'interval',
#'               intervalstart = 'lower')
#' plot.imputed(d, boxes = FALSE)
#' @export
plot.imputed <- function(x, boxes = FALSE, ...) {
  
  d <- x$proposals
  
  quantil <- ggplot2::ggplot(x$imputequantiles,
                    aes(x = seq(1:nrow(x$imputequantiles)),
                        ymin=q5,
                        lower=q25,
                        middle=q50,
                        upper=q75,
                        ymax=q95)) +
    geom_boxplot(stat="identity")+
    xlab('Observation')+
    ylab('Avg. quantiles of censored\n conditional bootmodel distribution')
  
  # Convert to Longformat
  d$observation <- 1:nrow(d)
  d <- reshape2::melt(d ,  id.vars = 'observation', variable.name = 'proposalVec')
  
  if (boxes) {
    imputations <- ggplot2::ggplot(data = d, aes(observation, value)) +
      geom_boxplot(aes(group = observation)) +
      xlab('Observation')+
      ylab('Proposals for observation [i]')
  }else {
    imputations <- ggplot2::ggplot() +
      geom_point(data = d, aes(observation, value)) +
      geom_point(data = data.frame(obs = 1:x$nreplacements,
                                   imputedx = x$imputedx), 
                 aes(x = obs, y = imputedx), color = 'red' )+
      xlab('Observation')+
      ylab('Proposals for observation [i]')
  }
  
  # Densities
  # Access censored covariate
  x1 <- as.character(x$mcall$xmu_formula[2])
  xobs <- x$Wobs$x1
  
  densities <-  ggplot() + 
    geom_density(data = data.frame(xobs), aes(x = xobs, color = "observed"), alpha = 0.4, size = 1.1) +
    geom_density(aes(x = value, y = ..density.., group = proposalVec, color = "proposal vector"), 
                 data = d, stat = "density", size = 1.1) +
    xlab('Covariate which includes defected data') +
    ylab('Density') +
    scale_color_discrete("")
  
  gridExtra::grid.arrange(quantil, imputations, densities, nrow = 1)
}

# description: y and defected not included.
# not a bug but a feature, that if only one var is remaining, its parallel lines

#' @title Andrews Curves of defected observations
#' @description Andrews Curves are a Fourier series upon the observations in
#'   data. They are a tool for detecting hidden groupings, and in this case of
#'   defected observations a tool for determining whether there is a clear
#'   structure in the remaining covariates, that may explain why a certain
#'   observation is likely to be defected. As it is an explorative tool, where
#'   the ordering of the variables determines the frequency that is affected
#'   respectively, it is highly recommended use various column orders. It may
#'   even be of use to some extent to employ Principle components. Note, that
#'   the defected, dependent and defect-indicator (and lowerbound in the
#'   intervalcase) variables are not considered for the andrews curve, as the
#'   information contained is ambigous and misleading. Particulary the dependent
#'   variable of the actual regression problem (not the imputation problem) is
#'   misleading, as it is caused by the covariates and not vice versa. Further
#'   Note, that given after deleting those columns only one covariate remains,
#'   the fourier will correctly return parallel lines: each value of that
#'   covariate is devided by sqrt(2). This is a feature not a bug.
#' @param dependent character. specifies the variable name of the dependent
#'   variable in the original regression problem (not the imputation problem)
#' @param ordering character vector, specifying the order of the variables in
#'   the andrews curve. Note that the ordering relates to the frequency in a
#'   fourier that is associated with a covariate.
#' @examples 
#' finterval = simulateData(n= 100,
#'                       param.formula = list(mu = ~exp(x1)+ x2+ x3, sigma = ~sin(x2)),
#'                       name = 'x1', subset = ~ (x2 < 0.3& x3<0.2), prob = 0.4, damage =list(c(0.3,0.9),c(1.2,1.5)),
#'                       family = 'NO',
#'                       correlation = matrix(c(1,0.3,0.2,
#'                                              0.3,1, 0.4,
#'                                              0.2,0.4,1), nrow = 3))
#' d = imputex(data = finterval$defected,
#'             xmu_formula= x1~ y+x2+x3,
#'             xsigma_formula = ~x2,
#'             xfamily = NO(mu.link = 'identity'),
#'             indicator = "indicator",
#'             censtype= 'interval',
#'             intervalstart = 'lower')
#' andrew(d, dependent = 'y')
#' @export
andrew <- function(object, dependent, ordering = NULL){
  
  if(class(object) != "imputed")
    stop("Argument 'object' has to be of class 'imputed'!")

  defected <- as.character(object$mcall$xmu_formula[[2]])
  data <- object$fulldata
  indicator <- object$mcall$indicator
  
  if(object$mcall$censtype == 'interval'){
    d <-  data[setdiff(names(data), c(defected, dependent, 'lower'))]
  }else{
    d <-  data[setdiff(names(data), c(defected, dependent))]
  }
  
  # reorder columns for later ease with positional matching in apply
  index <- which(names(d)== indicator)
  if(is.null(ordering)){
    d <- d[, c(setdiff(1:ncol(d), index),index)]
  } else {
    if(!identical(setdiff(names(d), c(ordering, indicator)) ,character(0))){
      stop('Names in Ordering were not found in provided data. All variable names must be specified')
    }
    d <- d[, c(ordering, indicator)]
  }
  
  if(length(setdiff(names(data), c(defected, dependent, indicator)))== 0){
    stop('dataframe must contain at least one variable apart from indicator, defected and dependent column')
  }
  
  #' @title Fourier series
  #' @param t axis position at which to evaluate
  #' @param obs observation vector
  #' @return Return of fourier series
  curveval <- function(t, obs){
    f <- obs[1] / sqrt(2)
    if(length(obs)>1){ 
      for (i in 2:length(obs)){
        if (i %% 2 == 0){ # even
          f <- f+ obs[i]*(sin((i-1)*t)+cos((i-1)*t))
        }else{ # odd
          f <- f+ obs[i]*(sin((i-2)*t)+cos((i-2)*t))
        }
      }
    }
    return(f)
  }
  
  p <- ggplot2::ggplot(data.frame(t = c(-pi, pi)), aes(t))
  p <- p + apply(
    d,
    MARGIN = 1,
    FUN = function(z)
      stat_function(
        fun = curveval,
        geom = "line",
        args = list(obs = z[1:(length(z)-1)]), # last is indicator
        color = z[length(z)] + 1  # based on indicator! color must be positive and dummy is 0/1
      )
  )
  print(p)
}

