#' @title Printing an object of class "imputed"
#' 
#' @param x Object of class "imputed".
#' @param ... print-specific arguments. See print() documentation.
#' 
#' @return Returns a print in the console.
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



#' @title Summarizing an object of class "imputed"
#' 
#' @param object Object of class "imputed".
#' @param ... summary-specific arguments. See summary() documentation.
#' 
#' @return Returns a summary in the console.
#' @export

summary.imputed <- function(object, ...) {
  
  # Access elements from result list in imputex:
  v <- object$nreplacements
  n <- object$nobservations
  sm <- summary(object$imputedx)
  mcall <- object$mcall
  cens_type <- object$censtype
  rounds <- object$m
  sr <- summary(object$distances)
  
  # Distance is undefined if interval or missing:
  r <- switch(is.null(object$distances) + 1, round(sr, 3), "Undefined")
  
  # Forming the shape of display:
  cat("\n", 
      cat("Call:  \n", paste(deparse(mcall), sep = "\n", collapse = "\n")),
      "\n",
      paste(round(v/n * 100, 2),"% of the observations are defected", sep = ""),
      "\n",
      paste("\n Number of observations:", n),
      paste("\n Type of censoring:", cens_type),
      paste("\n Number of proposals for each value to be replaced:", rounds),
      paste("\n Number of replacements:", v),
      "\n\n Imputed values:", 
      "\n",
      paste(names(sm), collapse = "     "),
      "\n",
      paste(round(sm, 3), collapse = "      "),
      "\n\n Distances of imputations to censorings:",
      "\n",
      paste(names(r), collapse = "     "),
      "\n",
      paste(r, collapse = "      ")
  )
}


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
#' @param boxes boolean. Indicating whether (2) should be displayed as a
#'   boxplot. Note that the median values are the imputation values.
#' @param ... plot-specific arguments. See plot() documentation.
#' 
#' @return Return figures with information on object of class "imputed".
#' @import ggplot2
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
  
  # Access censored covariate:
  x1 <- as.character(x$mcall$xmu_formula[2])
  xobs <- x$Wobs$x1
  
  densities <-  ggplot() + 
    geom_density(data = data.frame(xobs), 
                 aes(x = xobs, fill = "observed", color = "observed"),
                 alpha = 0.4, size = 0.5) +
    geom_density(data = d,
                 aes(x = value, 
                     y = ..density..,  
                     group = proposalVec, 
                     color = "proposal vector"), 
                 stat = "density", size = 0.5) +
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
#' @param ... Further arguments to be passed.
#'
#' @return Andrew's curves
#' @export

andrew <- function(object, ...) {
  UseMethod("andrew", object)
}


#' @title Andrew curve. Method of imputed 
#' @param object imputed.

andrew.imputed <- function (object, dependent, ordering = NULL){
  
  # read from object
  defected <- as.character(object$mcall$xmu_formula[[2]])
  data <- object$fulldata
  indicator <- object$mcall$indicator
  
  # corner case 'interval' has additional non informative column
  if(object$mcall$censtype == 'interval'){
    d <-  data[setdiff(names(data), c(defected, dependent, 'lower'))]
  }else{
    d <-  data[setdiff(names(data), c(defected, dependent))]
  }
  
  # reorder columns for andrew call/ and or via user shuffle
  index <- which(names(d)== indicator)
  if(is.null(ordering)){
    d <- d[, c(setdiff(1:ncol(d), index),index)]
  } else {
    if(!identical(setdiff(names(d), c(ordering, indicator)) ,character(0))){
      stop('Names in Ordering were not found in provided data. All variable names must be specified')
    }
    d <- d[, c(ordering, indicator)]
  }
  andrew(data = d)
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
#'   covariate is devided by sqrt(2). This is a feature not a bug.
#'
#' @param d dataframe. Contains observations rowwise and last column is a group
#'   indicator. The indicator is responsible for coloring of the curves.
#'   Notably, the input format exceeds the dummy format; any integer values can
#'   be used to indicate grouping.
#' @param t vector (sequence). At which the fourier series is to be evaluated.
#' 
#' @return ggplot of Andrews curve. Colored according to indicator.
#' @example
#' d = as.matrix(data.frame(x = c(1,2), y = c(3,4), ind = c(1,0)))
#' andrew(data = d)
#' d = as.matrix(data.frame(x = c(1,2,3), y = c(3,4,4), ind = c(1,0,2)))
#' andrew(data = d)

andrew <- function(data, t = seq(-pi, pi, length.out = 100)) {
  
  parameters <- data[, -ncol(data)]
  
  # fourier striped of parameters; a 'vector' (actually a list) filled with
  # unevaluated summands of fourier series, dependent on t without the parameter
  # factors. generates a 'vector' of unevaluated fourier series without the
  # parameter factors (observation's covariate values)
  # @param parameter matrix filled with fourier parameters of an observation rowwise
  stripedfourier <- function(parameters) {
    l = list(~1/sqrt(2))
    if(length(parameters)>1){
      for(i in 2:length(parameters)){
        if (i %% 2 == 0){ # even
          l[[i]] = as.formula(paste('~sin((' ,i, '-1)*t)+cos((', i, '-1)*t)'))
        }else{ # odd
          l[[i]] = as.formula(paste('~sin((',i,'-2)*t)+cos((',i, '-2)*t)' ))
        }
      }
    }
    return(l)
  }
  
  # executed only once: find the fourierseries of appropriate length.
  l <- stripedfourier(parameters[1,])
  
  # @description  Evaluate the Fourier series without the obseravational parameters
  # @param t vector. axis position(s) at which fourier series is to be evaluated
  # @param stripedfourier list. resulted fourier series expansion from stripedfourier()
  # @example
  # # d = as.matrix(data.frame(x = c(1,2), y = c(3,4), ind = c(1,0)))
  # # stripedt( t = c(1,2,3) ,  stripedfourier= l )
  stripedt = function(t, stripedfourier){
    
    # evaluate for all t (still list, as constant is evaluated only once and prevents
    # from coercing to matrix, as list is of unbalanced length)
    v = sapply(l, FUN= function(e) eval(e[[2]], envir = data.frame(t = t)))
    
    #! Workaround:
    
    # expand the constant to appropriate size
    v[[1]] = rep(v[[1]], length(t))
    
    # matrix of striped fourier (columnwise)
    v = sapply(v, function(x) unlist(x))
    
    # zip striped fourier with its parameter vector.
    return(v)
  }
  
  # executed only once to get the striped fouriervalues (dependent on t)
  v <- stripedt(t, stripedfourier= l)
  
  # @title zip fourier with the set of parameters
  # @description columnwise fully evaluated (all t values) fourier observations (actually only matrix product)
  # @param v stripedt matrix. The striped fourier series without considering observational
  # parameters. ('raw fourier')
  # @param parameters matrix. containing rowwise the covariates of each observation.
  zipfourier <- function(v, parameters){
    fourierobs = data.frame(v %*% t(parameters))
    return(fourierobs)
  }
  
  # convert to longformat
  observations <- zipfourier(v, parameters)
  observations$t <- t
  dat <- reshape2::melt(observations, id.vars = 't', variable.name = 'obs')
  
  # expand indicator of defected
  dat$indicator = rep(data[,ncol(data)], each= length(t))
  
  ggplot(data = dat, aes(x = t, y = value, group = obs, color = indicator))+
    geom_line()
}


