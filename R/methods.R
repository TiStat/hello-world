# Print S3 method -----------------------------------------------------------------------

#' @title Printing an object of class "imputed"
#' 
#' @param x Object of class "imputed".
#' 
#' @param ... print-specific arguments. See print() documentation.
#' 
#' @return Returns a print in the console.
#' 
#' @export

print.imputed <- function(x, ...) {
  
  # Access elements from result list in imputex:
  v <- x$nreplacements
  cens_type <- x$censtype
  n <- x$nobservations
  mcall <- x$mcall
  
  # Forming the shape of display:
  cat("\n", 
      cat("Call:"),
      paste(deparse(mcall), sep = "\n", collapse = "\n"),
      sprintf("\n\n Number of observations: %s\n", n),
      paste("Replaced", v, cens_type, "censored values") 
  )
}


# Summary S3 method ----------------------------------------------------------------------------

#' @title Summarizing an object of class "imputed"
#' 
#' @param object Object of class "imputed".
#' 
#' @param ... summary-specific arguments. See summary() documentation.
#' 
#' @return Returns a summary in the console.
#' 
#' @export

summary.imputed <- function(object, ...) {
  
  # Access elements from result list in imputex:
  v <- object$nreplacements
  n <- object$nobservations
  sm <- summary(object$imputedx)
  mcall <- object$mcall
  cens_type <- object$censtype
  
  # For missing or interval-censored data, distance is undefined:
  r <- ifelse(!is.null(object$distances), round(mean(object$distances), 3), "undefined")
  
  # Forming the shape of display:
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


# Plot S3 method --------------------------------------------------------------------------

#' @title Plot proposed imputations and their conditional censored
#'   distributions' average quantiles
#'   
#' @description Given an imputex object, this funciton plots: \cr (1) The
#'   approximate averaged quantiles of the valid part of the censored conditonal
#'   distributions, from which the proposed vectors were 'drawn' before the
#'   proposals were aggregated via median to become the imputed vector. \cr (2) The
#'   actual imputations for each observation. The red dots are the median
#'   proposals, which are in fact the final imputed vector. \cr (3) Density plots on
#'   each proposal vector. Consider that each proposal vector is drawn from a
#'   bootmodel-distribution, whose observations have their own parameters. So
#'   each proposalvector is a realization of a mixed distribution. Suppose, that
#'   the booted models are relatively similar, the density estimates on all
#'   proposal vectors (which are displayed) should be relatively similar with
#'   little variation.
#'   
#' @param x Object of class "imputex.
#' 
#' @param boxes boolean. Indicating whether (2) should be displayed as a
#'   boxplot. Note that the median values are the imputation values.
#'   
#' @param ... plot-specific arguments. See plot() documentation.
#' 
#' @examples 
#' # Simulating a dataset
#' rinterval = simulateData(n= 300,
#'                          param.formula = list(mu = ~exp(x1), sigma = ~sin(x2)),
#'                          name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 ,
#'                          damage =list(c(0.8, 0.99), c(1.2,1.5)),
#'                          family = 'NO',
#'                          correlation = NULL)
#'                          
#' d  <- imputex(data = rinterval$defected,
#'               xmu_formula= x1~y,
#'               indicator = "indicator",
#'               censtype = 'interval',
#'               intervalstart = 'lower')
#'               
#' plot(d, boxes = FALSE)
#' 
#' @return Return figures with information on object of class "imputed".
#' 
#' @export

plot.imputed <- function(x, boxes = FALSE, ...) {
  
  # Access elements from result list in imputex:
  d <- x$proposals
  im <- x$imputequantiles
  
  quantil <- ggplot(im, aes(x = seq(1:nrow(im)),
                            ymin=q05,
                            lower=q25,
                            middle=q50,
                            upper=q75,
                            ymax=q95)) +
    
    geom_boxplot(stat="identity") +
    
    xlab('Observation') +
    ylab('Avg. quantiles of censored\n conditional bootmodel distribution') +
    theme(axis.title=element_text(size=11,face="bold"))
  
  # Convert to Longformat:
  d$observation <- 1:nrow(d)
  d <- reshape2::melt(d ,  id.vars = 'observation', variable.name = 'proposalVec')
  
  if (boxes == TRUE) {
    
    imputations <- ggplot(data = d, aes(observation, value)) +
      geom_boxplot(aes(group = observation)) +
      
      stat_summary(geom = "crossbar", width = 0.65, color="red",
                   fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
      
      xlab('Observation') +
      ylab('Proposals for observation [i] & median') +
      theme(axis.title=element_text(size=11,face="bold"))
    
  } else {
    imputations <- ggplot() +
      geom_point(data = d, aes(observation, value)) +
      
      geom_point(data = data.frame(obs = 1:x$nreplacements, imputedx = x$imputedx),
                 aes(x = obs, y = imputedx), color = 'red') +
      
      xlab('Observation') +
      ylab('Proposals for observation [i] & median') +
      theme(axis.title=element_text(size=11,face="bold"))
  }
  
  # Densities;
  # Access censored covariate:
  x1 <- as.character(x$mcall$xmu_formula[2])
  xobs <- x$Wobs$x1
  
  densities <-  ggplot() + 
    geom_density(data = data.frame(xobs), 
                 aes(x = xobs, fill = "observed", color = "observed"),
                 alpha = 0.4, size = 0.5) +
    
    geom_density(aes(x = value, y = ..density..,  
                     group = proposalVec, 
                     color = "proposal vector"), 
                 data = d, stat = "density", size = 0.5) +
    
    scale_fill_discrete(guide=FALSE) +
    xlab('Covariate which includes defected data') +
    ylab('Density') +
    scale_color_discrete("") +
    theme(axis.title=element_text(size=11,face="bold"))
  
  # Display plots in one window:
  gridExtra::grid.arrange(quantil, imputations, densities, nrow = 1)
}








# Andrew's curves ----------------------------------------------------------------------------



#' @title Andrew's curves
#' 
#' @description `andrew` is a generic function used to display Andrew's curves of the
#' independent variables of an object of class "imputed". The function has just one single
#' method. See documentation of `andrew.imputed` for further description.
#'
#' @param object Object of class "imputed".
#' 
#' @param ... Further arguments to be passed.
#'
#' @return Andrew's curves
#' 
#' @export

andrew <- function(object, ...) {
  UseMethod("andrew", object)
}





#' @title Andrew's Curves of defected observations
#' 
#' @description Andrew's Curves are a Fourier series upon the observations in
#'   data. They are a tool for detecting hidden groupings, and in this case of
#'   defected observations, a tool for determining whether there is a clear
#'   structure in the remaining covariates, that may explain why a certain
#'   observation is likely to be defected. As it is an explorative tool, where
#'   the ordering of the variables determines the frequency that is affected
#'   respectively, it is highly recommended to use various column orders. It may
#'   even be of use to some extent to employ Principle Components. Note, that
#'   the defected, dependent and defect-indicator (and lower bound in the
#'   interval case) variables are not considered for the Andrew's curve, as the
#'   information contained is ambigous and misleading. Particulary, the dependent
#'   variable of the actual regression problem (not the imputation problem) is
#'   misleading, as it is caused by the covariates and not vice versa. Further,
#'   note that after deleting those columns only one covariate remains, so
#'   the fourier will correctly return parallel lines: each value of that
#'   covariate is devided by sqrt(2). This is a feature, not a bug.
#'   
#' @param object Object of class "imputed".
#' 
#' @param dependent character. Specifies the variable name of the dependent
#'   variable in the original regression problem (not the imputation problem).
#'   
#' @param ordering character vector, specifying the order of the variables in
#'   the Andrew's curve. Note that the ordering relates to the frequency in a
#'   fourier that is associated with a covariate.
#'   
#' @param ... Further arguments to be passed.
#'   
#' @return Returns Andrew's curves figure.
#'   
#' @examples 
#' finterval = simulateData(n= 100,
#'                          param.formula = list(mu = ~exp(x1) + x2 + x3, 
#'                          sigma = ~sin(x2)),
#'                          name = 'x1', 
#'                          subset = ~ (x2 < 0.3 & x3 < 0.2),
#'                          prob = 0.4, 
#'                          damage =list(c(0.3, 0.9), c(1.2, 1.5)),
#'                          family = 'NO',
#'                          correlation = matrix(c(1, 0.3, 0.2, 0.3, 1, 0.4, 0.2, 0.4, 1), nrow = 3))
#' 
#' d <- imputex(data = finterval$defected,
#'              xmu_formula= x1~ y + x2 + x3,
#'              xsigma_formula = ~x2,
#'              xfamily = NO(mu.link = 'identity'),
#'              indicator = "indicator",
#'              censtype= 'interval',
#'              intervalstart = 'lower')
#'             
#' andrew(d, dependent = 'y')
#' 
#' @export

andrew.imputed <- function(object, dependent, ordering = NULL, ...) {
  
  if(class(object) != "imputed")
    stop("Argument 'object' has to be of class 'imputed'!")
  
  # Access some elements from imputex:
  defected <- as.character(object$mcall$xmu_formula[[2]])
  data <- object$fulldata
  indicator <- object$mcall$indicator
  
  if(object$mcall$censtype == 'interval'){
    d <-  data[setdiff(names(data), c(defected, dependent, 'lower'))]
  }else{
    d <-  data[setdiff(names(data), c(defected, dependent))]
  }
  
  # Reorder columns for later ease with positional matching in apply:
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
    stop('data frame must contain at least one variable apart from indicator, defected and dependent column')
  }
  
  # Following function evaluates fourier series at "t" for parameter set "obs", where
  # t is the axis position at which to evaluate and obs is the observation vector:
  curveval <- function(t, obs){
    f <- obs[1] / sqrt(2)
    if(length(obs)>1){ 
      for (i in 2:length(obs)){
        if (i %% 2 == 0){ # even
          f <- f + obs[i]*(sin((i-1)*t)+cos((i-1)*t))
        }else{ # odd
          f <- f + obs[i]*(sin((i-2)*t)+cos((i-2)*t))
        }
      }
    }
    return(f)
  }
  
  p <- ggplot(data.frame(t = c(-pi, pi)), aes(t))

  p <- p + apply(d, MARGIN = 1, FUN = function(z) stat_function(fun = curveval,
                                                                geom = "line",
                                                                args = list(obs = z[1:(length(z)-1)]), # last is indicator
                                                                # based on indicator! color must be positive and dummy is 0/1
                                                                color = z[length(z)] + 1)
  )
    
  
  print(p)
}

