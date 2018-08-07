#' Title
#'
#' @param ymu_formula
#' @param ysigma_formula
#' @param ynu_formula
#' @param ytau_formula
#' @param yfamily
#' @param xmu_formula
#' @param xsigma_formula
#' @param xnu_formula
#' @param xtau_formula
#' @param xfamily
#' @param data
#' @param indicator
#' @param censtype
#' @param ...
#'
#' @return
fitfull = function(ymu_formula,
                   ysigma_formula,
                   ynu_formula,
                   ytau_formula,
                   yfamily,
                   xmu_formula,
                   xsigma_formula, # specify default here?
                   xnu_formula,
                   xtau_formula,
                   xfamily,
                   data,
                   indicator,
                   censtype,
                   ...){      # seperate ellipsis for y and for x fit!!
  imputation = imputex(xmu_formula,
                       xsigma_formula = ~1,
                       xnu_formula = ~1,
                       xtau_formula = ~1,
                       xfamily = NO(mu.link = 'identity'),
                       data,
                       indicator,
                       censtype,
                       ...)

  yfit = gamlss(formual = ymu_formula,
                sigma.formula = sigma_formula,
                nu.formula = ynu_formula,
                tau.formula = ytau_formula,
                family = yfamily,
                ydata = imputation$fulldata )
}



