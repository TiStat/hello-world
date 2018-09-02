#' @title Andrew's curves
#' 
#' @description `andrew` is a generic function used to display Andrew's curves of the
#' independent variables of an object of class "imputed". The function has just one single
#' method. See documentation of `andrew.imputed` for further description.
#'
#' @param object Object of class "imputed".
#' @param ... Further arguments to be passed to methods.
#'
#' @return Andrew's curves
#' @export

andrew <- function(object, ...) {
  UseMethod("andrew", object)
}


#' @title Andrew's Curves of covariates for "imputed" class objects
#' 
#' @description Andrew's Curves are a Fourier series upon the observations in
#'   data. They are a tool for detecting hidden groupings and in this case of
#'   defected observations, they're a tool for determining whether there is a clear
#'   structure in the remaining covariates, that may explain why a certain
#'   observation is likely to be defected. As it is an explorative tool, where
#'   the ordering of the variables determines the frequency that is associated
#'   respectively, it is highly recommended to use various column orders. It may
#'   even be of use to some extent to employ Principle Components. Note that
#'   the defected, dependent and defect-indicator (and lower bound in the
#'   interval case) variables are not considered for the Andrew's curve, as the
#'   information contained is ambiguous and misleading. Particulary, the dependent
#'   variable of the actual regression problem (not the imputation problem) is
#'   misleading, as it is caused by the covariates and not vice versa. Further,
#'   note that after deleting those columns only one covariate remains, so
#'   the fourier will correctly return parallel lines: each value of that
#'   covariate is devided by sqrt(2). This is a fourier feature, not a bug.   
#'   
#' @param object Object of class "imputed".
#' @param dependent character. Name of the dependent variable in the original
#'   (not imputation) regression problem. It is removed, as the information
#'   contained is dubious: Covariates cause the dependent and not vice versa.
#' @param ordering vector of characters. Names of the covariates supplied to
#'   imputex. The argument is optional and allows to shuffle the dataframe.
#'   Thereby, the covariates are associated with different Fourier frequencies.
#'   It is highly recommended to make use of this option. As syntax sugar, it is
#'   possible to specify only the first few variables and leave the remaining
#'   ordering in the dataframe intact.
#' @param ... Further arguments to be passed (e.g. of the imputed object)
#'   
#' @examples 
#' fright <- simulateData(n= 150,
#'                       param.formula = list(mu = ~ exp(x1) + x2+ x3, sigma = ~sin(x2)),
#'                       name = 'x1', subset = ~ x1 > 0.6, prob = 0.8 , damage = 1/3,
#'                       family = 'NO',
#'                       correlation = matrix(c(1, 0.3, 0.2,
#'                                              0.3, 1, 0.4,
#'                                              0.2, 0.4, 1), nrow = 3))
#'                                             
#' d <- imputex(data = fright$defected,
#'              xmu_formula = x1 ~ y + x2 + x3,
#'              xsigma_formula = ~x2,
#'              xnu_formula = ~1,
#'              xtau_formula = ~1,
#'              xfamily = NO(mu.link = 'identity'),
#'              indicator = "indicator",
#'              censtype= 'right' )
#'              
#' andrew(object = d, dependent = 'y', ordering = c('x3'))
#' @export


andrew.imputed <- function (object, dependent, ordering = NULL, ...) {
  
  # Read from object:
  defected <- as.character(object$mcall$xmu_formula[[2]])
  data <- object$fulldata
  indicator <- object$mcall$indicator
  
  # Corner case 'interval' has additional non informative column:
  if (object$mcall$censtype == 'interval') {
    d <- data[setdiff(names(data), c(defected, dependent, 'lower'))]
  } else{
    # General case, removes defected & dependent
    d <- data[setdiff(names(data), c(defected, dependent))]
  }
  
  # Reorder/Shuffle columns
  # by default, indicator column is set as last.
  shuffle <- function(d, ordering, indicator) {
    index <- which(names(d) == indicator)
    if (is.null(ordering)) {
      d <- d[ , c(setdiff(1:ncol(d), index), index)]
      
    } else if (!identical(setdiff(names(d), c(ordering, indicator)) , character(0))) {
      remainder <- setdiff(names(d), c(ordering, indicator))
      d <- d[ , c(ordering, remainder, indicator)]
      
    } else {
      d <- d[ , c(ordering, indicator)]
    }
  }
  
  d <- shuffle(d, ordering, indicator)
  
  # Function 'andrewcore' scripted below.
  andrewcore(data = d)
}

#' @title Core function to be passed to andrew.imputed
#'
#' @description `andrew.imputed` prepares the data frame to be passed to this function.
#'
#' @param data data.frame (or matrix). Contains observations rowwise and last column is a
#'   group indicator. The indicator is responsible for coloring of the curves.
#'   Notably, the input format exceeds the dummy format; any integer values can
#'   be used to indicate grouping.
#' @param t vector (sequence). At which the fourier series is to be evaluated.
#' 
#' @return Plot of Andrew's curve. Colored according to indicator.

andrewcore <- function(data, t = seq(-pi, pi, length.out = 100)) {
  
  # reduce dataframe to parameter part & and prevent R from coercing to vector
  # if only one covariate remains after removal of indicator
  parameters <- data[ , -ncol(data), drop = FALSE]
  
  # fourier striped of parameters; list (as surrogate for a functional vector) filled with
  # unevaluated summands of fourier series, dependent on t without the parameter
  # factors. 
  # @param nparameter Number of Parameters of dataframe, for which the fourier is expanded.
  # @note The Workaround t/t is due, as eval of v would otherwise not expand
  #   the constant to appropriate length and return an unbalanced list, which in
  #   turn will not be unlisted in a matrix
  stripedfourier <- function(nparameters) {
    l = list(~ 1/sqrt(2) * t/t) # unfortunate workaround *1
    if(nparameters > 1) {
      for(i in 2:nparameters) {
        if (i %% 2 == 0){ # even
          l[[i]] <- as.formula(paste('~sin((' ,i, '-1)*t) + cos((', i, '-1)*t)'))
        }else{ # odd
          l[[i]] <- as.formula(paste('~sin((',i,'-2)*t) + cos((',i, '-2)*t)' ))
        }
      }
    }
    return(l)
  }
  
  # Fourier Series' summands striped of observational parameter of appropriate
  # length. unevaluated and dependent on t
  l <- stripedfourier(nparameters = length(parameters[1, ]))
  
  # Evaluate the Fourier series' summands without the obseravational parameters for t 
  # result is striped (of parameters) Fourier matrix. (summands evaluated at t rowwise)
  v <- sapply(l, FUN= function(e) eval(e[[2]], envir = data.frame(t = t)))
  
  # zip fourier: scale the raw summands evaluated at t with all of their
  # respective observational parameter vectors. result is columnwise the fourier
  # expansion for each row of the original data. Note that t() implicitly
  # coerces to matrix and allows argument 'parameters' in andrew() to be either
  # matrix or dataframe. Coerce to Dataframe for ggplot.
  fourierobs <- data.frame(v %*% t(parameters))
  
  # convert to longformat to arrive at an automated color scheme
  fourierobs$t <- t
  dat <- reshape2::melt(fourierobs, id.vars = 't', variable.name = 'obs')
  
  # expand indicator of defected
  dat$indicator <- rep(data[ , ncol(data)], each= length(t))
  
  ggplot(data = dat, aes(x = t, y = value, group = obs, color = factor(indicator))) +
    geom_line() +
    scale_color_manual(name = "Observation", 
                       labels = c("Fully observed", "Defected"),
                       values = c("black", "red"))
}


