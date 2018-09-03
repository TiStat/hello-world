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
  impvar <- summary(object$imputevariance)
  
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
      paste("\n Number of proposals for each defected observation:", rounds),
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
      paste(r, collapse = "      "),
      "\n\n Imputation variances:",
      "\n",
      paste(names(impvar), collapse = "     "),
      "\n",
      paste(round(impvar, 3), collapse = "      ")
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





