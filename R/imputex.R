# Dependencies -----------------------------------------------------------------
# library('ggplot2')
library(gamlss)

#'@description split dataset
#'@param data dataframe.
#'@param indicator character. indicating the name of the missing/ (right/left)
#' censored observation dummy variable in data
W <- function(data, indicator){ # function name als handlungsanweisung formulieren
  df_obs <- data[data[indicator] == 0, ]
  df_cens <- data[data[indicator] == 1, ]
  return(list(obs = df_obs, cens = df_cens))
}



# Algorithm --------------------------------------------------------------------
#' Title
#'
#' @Description
#'
#' @param data data.frame containing a dummy censoring indicator, 0 if not
#'   indicator, 1 if indicator
#' @param indicator character. Name of dummy column in data, which indicates the
#'   damaged observation.
#' @param censtype character. type of the damaged observation 'missing', 'right', 'left', 'interval'
#' @param xmu_formula
#' @param xsigma_formula
#' @param xnu_formula
#' @param xtau_formula
#' @param xfamily gamlss family object.
#' @param ... additional arguments passed in the respective gamlss fit # check if these arguments are available in gamlss

#'
#' @return
imputex <- function(xmu_formula,
                    xsigma_formula = ~1,
                    xnu_formula = ~1,
                    xtau_formula = ~1,
                    xfamily = NO(mu.link = 'identity'),
                    data,
                    indicator,
                    censtype,
                    ...)
{
  if(!is.data.frame(data)){ #| is.empty(d)
    stop('data must be data.frame')
  }

  if(!((is.character(indicator)) & (indicator %in% names(data)))){
    stop('indicator must be a column name in data')
  }

  # save for user call in summary
  mcall <- match.call()

  # split dataset in fully observed & missing/censored data
  Wdat <- W(data, indicator)
  censor = as.character(xmu_formula[[2]])

  # Step 1: fit gamlss with user specified xfamily and formula on observed data
  obsmodel <- gamlss(
    formula = xmu_formula,
    sigma_formula = xsigma_formula,
    nu_formula = xnu_formula,
    tau_formula = xtau_formula,
    family = xfamily,
    data = Wdat$obs,
    ...)

  # Step 2: Resampling from fitted model
  # note that these are independent draws from the same distribution. Reframe to m vectors:
  draws <- family_fun(object = obsmodel,
                      func = "r",
                      predictdata = Wdat$obs,
                      n = nrow(Wdat$obs)*nrow(Wdat$cens),
                      fitdata = Wdat$obs)
  draws <- data.frame(matrix(draws,
                             nrow = nrow(Wdat$obs),
                             ncol = nrow(Wdat$cens)))
  # Basis for respective Bootstap samples (draws[,j], W).
  drawsW <- cbind(draws, Wdat$obs)

  # Bootstap samples on simulated vectors (rowwise)
  boot <- drawsW[sample(x= 1:nrow(drawsW),
                        size = nrow(Wdat$obs),
                        replace = TRUE), ]

  # step 3 estimate param. on each respective set {x*boot(j), W_obs} for all j
  bootmodel <- list()
  imputemat <- data.frame(X1 = vector(length= nrow(Wdat$cens)))
  for (i in 1:ncol(draws)){
    # iterate only over the names of booted vectors.

    # manipulate the formula for each estimation to get the respective
    # booted x column regressed on W
    bootformula <- as.formula(paste(names(boot)[i], '~',
                                    as.character(as.vector(xmu_formula)[3]),
                                    sep = ''))
    bootmodel[[i]] <- gamlss(
      formula = bootformula,
      sigma_formula = xsigma_formula,
      nu_formula = xnu_formula,
      xtau_formula = xtau_formula,
      family = xfamily,
      data = boot,...)

    # Simulate data from the corresponding fitted distribution.

    imputecandidate <- names(boot)[i]
    imputemat[[imputecandidate]] = samplecensored(
      object = bootmodel[[i]],
      censtype,
      predictdata = Wdat$cens,
      censor,
      fitdata = boot
    )
  }
  imputemat$imputedx <- apply(imputemat, MARGIN = 1, mean)

  # complete data with imputations
  Wdat$cens[censor] = imputemat$imputedx
  fulldata = rbind(Wdat$obs,Wdat$cens)

  return(list(imputations = imputemat, fulldata = fulldata)) # edit output format!
}

d2 <- imputex(data= data,
             xmu_formula= x1~y+x2,
             xsigma_formula = ~1,
             xnu_formula = ~1,
             xtau_formula = ~1,
             xfamily = NO(mu.link = 'identity'),
             indicator = "indicator",
             censtype = 'right')

# Mit ellipsis geändert
d3 <- imputex(data= data,
             xmu_formula= x1~y+x2,
             xsigma_formula = ~1,
             xnu_formula = ~1,
             xtau_formula = ~1,
             xfamily = NO(mu.link = 'identity'),
             indicator = "indicator",
             censtype = 'right',method = CG())
