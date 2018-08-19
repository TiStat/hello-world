# Create print method for "imputed" class
#' Title
#'
#' @param x 
#'
#' @return What should be printed for "imputed" class object
#' @export
print.imputed <- function(x){
  
  mcall <- noquote(capture.output(x$mcall)) #call
  m <- x$number_of_imputations
  cens_type <- x$censoring_type
  n <- x$number_of_observations
  
  cat("\n", 
      cat("Call:"),
      mcall,
      sprintf("\n\n Number of observations: %s\n", n),
      paste("Imputed", m, cens_type, "censored values") )
}


#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples 
summary.imputed <- function(x) {
  
  m <- x$number_of_imputations #number of imputations
  n <- x$number_of_observations
  #  v <- x$variance #discrepancy
  #  aic <- x$AIC
  mcall <- noquote(capture.output(x$mcall)) #call
  cens_type <- x$censoring_type
  
  cat("\n", 
      cat("Call:  \n", mcall),
      paste("\n Number of observations:", n),
      paste("\n Type of censoring:", cens_type),
      paste("\n Number of imputations:", m))
}

#' @title plot imputations and their quantiles
#' @description Given an imputex object, this funciton plots the imputations and
#'  optionally plots the approximate averaged quantiles of the valid part of the
#'  censored conditonal distributions, from which the proposed vectors were
#'  'drawn', before the proposals were aggregated to become the imputed vector.
#' @param object imputex object.
#' @example plotimputations(d, boxes = FALSE, quantiles = TRUE)
plot.imputed <- function(object, boxes = TRUE, quantiles = FALSE) {

  d = object$imputations
  
  if(quantiles == TRUE){
    plot1 <- ggplot(object$impquantiles,
                    aes(x = seq(1:nrow(object$impquantiles)),
                        ymin=q5,
                        lower=q25,
                        middle=q50,
                        upper=q75,
                        ymax=q95)) +
      geom_boxplot(stat="identity")+
      xlab('Observation')+
      ylab('Avg. quantiles of censored\n conditional bootmodel distribution')
    

  }
  
  
  # Convert to Longformat
  d$observation = seq(1, nrow(d))
  d <- melt(d ,  id.vars = 'observation', variable.name = 'proposalVec')
  
  if (boxes) {

    plot2 <- ggplot(data = d, aes(observation, value)) +
      geom_boxplot(aes(group = observation)) +
      ylab('Proposals for observation [i]') # NOTE that red dots are based on mean, boxes display median.
  }else {
    plot2 <- ggplot() +
      geom_point(data = d, aes(observation, value)) +
      geom_point(data = data.frame(obs = 1:length(object$imputedx),imputedx = object$imputedx), aes(x = obs, y = imputedx), color = 'red' )+
      ylab('Proposals for observation [i]')
  }
  
  if (exists("plot1")) {
    grid.arrange(plot1, plot2, ncol=2) # make yaxis equal!

  } else {
    plot2
  }
}

#' @param data data.frame with defected observations 
#' @param defected character. specifys which column is defected.
#' @param indicator character. giving the name of the dummy variable, 
andrew = function(data, defected, indicator){
  
  if(length(names(data) != indicator)<=2){
    stop('dataframe must contain at least two variables apart from indicator and defected column')
  }
  
  # reorder columns for later ease with positional matching in apply
  d =  data[names(data) != defected]
  index = which(names(d)== indicator)
  d = d[, c(setdiff(1:ncol(d), index),index)]
  
  #' @param obs observation vector
  curveval = function(t, obs){
    f = obs[1] / sqrt(2)
    if(length(obs)>1){ # if branch will never be jumped over (due to stop!), but correct fourier function
      for (i in 2:length(obs)){
        if (i %% 2 == 0){ # even
          f = f+ obs[i]*(sin((i-1)*t)+cos((i-1)*t))
        }else{ # odd
          f = f+ obs[i]*(sin((i-2)*t)+cos((i-2)*t))
        }
      }
    }
    return(f)
  }
  
  p = ggplot(data.frame(t = c(-pi, pi)), aes(t))
  p = p + apply(
    d,
    MARGIN = 1,
    FUN = function(z)
      stat_function(
        fun = curveval,
        geom = "line",
        args = list(obs = z[1:length(z)]),
        color = z[length(z)] + 1  # must be positive and dummy is 0/1
      )
  )
  print(p)
}

