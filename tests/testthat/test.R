# Define some values from simulateData.R for unit test realisation --------------------------------------------

# Right censored data
rd <- simulateData(n= 300,
                   param.formula = list(mu = ~exp(x1) + x2 + x3, sigma = ~sin(x2)),
                   name = 'x1', subset = ~ (x2 > 0.6 & x3 > 0.7), prob = 0.75, 
                   damage =c(0.3, 0.9), family = 'NO', 
                   correlation = NULL)$defected
# Left censored data
ld <- simulateData(n= 300,
                   param.formula = list(mu = ~exp(x1) + x2 + x3, sigma = ~sin(x2)),
                   name = 'x1', subset = ~ (x2 < 0.3 & x3 < 0.4), prob = 0.8, 
                   damage =c(0.3, 0.9), family = 'NO', 
                   correlation = NULL)$defected
# Interval censored data
id <- simulateData(n= 300, 
                   param.formula = list(mu = ~exp(x1) + x2 + x3, sigma = ~sin(x2)),
                   name = 'x1', subset = ~ (x2 < 0.3 & x3 < 0.2), prob = 0.4, 
                   damage =list(c(0.3, 0.9), c(1.2, 1.5)), family = 'NO',
                   correlation = NULL)$defected

# dot(.) - indicator to prevent multicollinearity!
rmodel <- gamlss(formula = y ~ . -indicator, data=rd) 
lmodel <- gamlss(formula = y ~ . -indicator, data=ld)
imodel <- gamlss(formula = y ~ . -indicator, data=id)

# Number of censored values. This is important for the unit tests!
nr <- length(rd$x1[rd$indicator==1]) 
nl <- length(ld$x1[ld$indicator==1])
ni <- length(id$x1[ld$indicator==1])

#data frame to predict on, not a prediction!!!
rpredict.df <- data.frame(x1 = runif(n = nr), x2 = runif(n = nr), x3 = runif(n = nr), indicator = 1) 
lpredict.df <- data.frame(x1 = runif(n = nl), x2 = runif(n = nl), x3 = runif(n = nl), indicator = 1)
ipredict.df <- data.frame(x1 = runif(n = ni), x2 = runif(n = ni), x3 = runif(n = ni), indicator = 1)

rimpute <-  imputex(xmu_formula = x1 ~ y + x2 + x3, data = rd, indicator = "indicator", censtype = "right")
limpute <-  imputex(xmu_formula = x1 ~ y + x2 + x3, data = ld, indicator = "indicator", censtype = "left")
iimpute <-  imputex(xmu_formula = x1 ~ y + x2 + x3, data = id, indicator = "indicator", censtype = "interval", intervalstart = "lower")


# Unit tests for family_fun -----------------------------------------------------------------------------------
context('Evaluate family functions')

test_that('Test that n is a multiple of nrow(predictdata)', {
  expect_error(family_fun(object = rmodel, fitdata = rd, predictdata = rpredict.df, func = 'r', n = nrow(rpredict.df) + 1))
})

test_that('Test that mismatching arguments are stoped',{
  expect_error(family_fun(rmodel, func = 'r',rd, rpredict.df, x = 0.5))
})

test_that('Test if output is as expected', {
  expect_is(family_fun(rmodel, func = 'r',rd, rpredict.df, n = nrow(rpredict.df)), 'numeric')
  expect_is(family_fun(rmodel, func = 'p',rd, rpredict.df, q = c(0.5, 0.7)), 'numeric')
  expect_equal(length(family_fun(rmodel, func = 'r',rd, rpredict.df, n = nrow(rpredict.df))), nrow(rpredict.df))
})

test_that("Test that there are no NA's for censored data",{
  expect_equal(nrow(rpredict.df),length(family_fun(rmodel, func = 'd', rd, rpredict.df, x = 5)))
  # Next is only valid, if func is not 'r', since they are called at different seeds!
  expect_identical(family_fun(rmodel, func = 'd', rd, rpredict.df, x = 10), na.omit(family_fun(rmodel, func = 'd', rd, rpredict.df, x = 10)))
})


# Unit tests for samplecensored ----------------------------------------------------------------------------
context('Evaluate inverse sampling in samplecensored')

if (nrow(rpredict.df) == 1) {
  test_that('Test that quantiles of samplecensored are of the right class', {
    expect_is(samplecensored(rmodel ,censtype = 'right', rpredict.df, rd, censor = "x1")$quantiles, 'numeric')
  })
  
} else {
  test_that('Test that quantiles of samplecensored are of the right class', {
    expect_is(samplecensored(rmodel ,censtype = 'right', rpredict.df, rd, censor = "x1")$quantiles, 'matrix')
  })
}


# Unit tests for imputex ---------------------------------------------------------------------------------
context('Some tests for imputex')

test_that('Test that return is as expected', {
  expect_is(rimpute$proposals , 'data.frame') 
  expect_identical(rimpute$proposals, na.omit(rimpute$proposals))
  expect_equal(ncol(rimpute$proposals), rimpute$m)
  expect_is(rimpute$fulldata , 'data.frame')
  expect_equal(nrow(rimpute$fulldata), nrow(rd))
  expect_equal(ncol(rimpute$fulldata), ncol(rd))
  expect_is(rimpute$imputequantiles , 'data.frame')
  expect_equal(nrow(rimpute$imputequantiles), rimpute$nreplacements)
  expect_identical(rimpute$imputequantiles, na.omit(rimpute$imputequantiles))
  expect_true(all(rimpute$imputevariance > 0))
})

test_that('Test that call stops correctly with the valid error message', {
  rd$indicator = NULL
  expect_error(imputex(xmu_formula = x1 ~ y + x2, data = rd, indicator = "indicator", censtype = "right"), 'indicator must be a column name in data')
})

# right censored values and imputations
rcens_values <- rd$x1[rd$indicator == 1] # censored values
rimps <- rimpute$imputedx     # imputed values

# left censored values and imputations
lcens_values <- ld$x1[ld$indicator == 1] # censored values
limps <- limpute$imputedx     # imputed values

# interval censored values and imputations
iupper <- id$x1[id$indicator == 1] # upper bound censoring
ilower <- id$lower[id$indicator == 1] # lower bound censoring
iimps <- iimpute$imputedx     # imputed values

test_that('Valid draws, i.e. they are at least drawn from the valid region',{
  expect_true(all(rimps >= rcens_values)) # x1 defines lower bound!
  expect_true(all(limps <= lcens_values)) # x1 defines upper bound!
  expect_true(all(iimps >= ilower & iimps <= iupper)) # x1 defines upper bound and "lower" is the lower bound!
})

