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
  df = object$imputations
  
  if(quantiles == TRUE){
   plot1 <- ggplot(object$impquantiles,
                 aes(x = seq(1:nrow(object$impquantiles)),
                     ymin=q5,
                     lower=q25,
                     middle=q50,
                     upper=q75,
                     ymax=q95)) +
            geom_boxplot(stat="identity")+
            xlab('draw')+
            ylab('quantiles')

  }
  
  
  # Convert to Longformat
  df$observation = seq(1, nrow(df))
  df <- melt(df ,  id.vars = 'observation', variable.name = 'proposalVec')
  
  if (boxes) {
    plot2 <- ggplot(data = subset(df, df$proposalVec != 'imputedx'), aes(observation, value)) +
             geom_boxplot(aes(group = observation)) +
             geom_point(data = subset(df, df$proposalVec == 'imputedx'), aes(observation, value, color = 'red'))+
             ylab('Proposals for observation [i]') # NOTE that red dots are based on mean, boxes display median.
  }else {
   plot2 <- ggplot() +
             geom_point(data = subset(df, df$proposalVec != 'imputedx'), aes(observation, value)) +
             geom_point(data = subset(df, df$proposalVec == 'imputedx'), aes(observation, value, color = proposalVec))+
             ylab('Proposals for observation [i]')
  }
  
  if (exists("plot1")) {
  grid.arrange(plot1, plot2, ncol=2)
  } else {
    plot2
  }
}


