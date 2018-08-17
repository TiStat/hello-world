# Define some values from simulateData.R for unit test realisation --------------------------------------------
rd <- simulateData(n = 100)$defected
model <- gamlss(formula = y ~ ., data=rd)
nl <- length(rd$x1[rd$indicator==1]) # number of censored values. This is important for the unit tests!
predict.df <- data.frame(x1 = runif(n = nl), x2 = runif(n = nl), indicator = 1) #data frame to predict on, not a prediction!!!
impute <-  imputex(xmu_formula = x1 ~ y + x2, data = rd, indicator = "indicator", censtype = "right")

# Unit tests for family_fun -----------------------------------------------------------------------------------
context('Evaluate family functions')

test_that('Test that mismatching arguments are stoped',{
  expect_error(family_fun(model, func = 'r',rd, predict.df, x = 0.5))
})

test_that('Test if output is as expected', {
  expect_is(family_fun(model, func = 'r',rd, predict.df, n = 10), 'numeric')
  expect_is(family_fun(model, func = 'p',rd, predict.df, q = c(0.5, 0.7)), 'numeric')
  expect_equal(length(family_fun(model, func = 'r',rd, predict.df, n = 10)), 10)
})

test_that("Test that there are no NA's for censored data",{
  expect_equal(nrow(predict.df),length(family_fun(model, func = 'd', rd, predict.df, x = 5)))
  # Next is only valid, if func is not 'r', since they are called at different seeds!
  expect_identical(family_fun(model, func = 'd', rd, predict.df, x = 10), na.omit(family_fun(model, func = 'd', rd, predict.df, x = 10)))
})

# Unit tests for samplecensored ----------------------------------------------------------------------------
context('Evaluate inverse sampling in samplecensored')

if (nrow(predict.df) == 1) {
  test_that('Test that quantiles of samplecensored are of the right class', {
    expect_is(samplecensored(model ,censtype = 'right', predict.df, rd, censor = "x1")$quantiles, 'numeric')
  })
  
} else {
  test_that('Test that quantiles of samplecensored are of the right class', {
    expect_is(samplecensored(model ,censtype = 'right', predict.df, rd, censor = "x1")$quantiles, 'matrix')
  })
}

# Unit tests for imputex ---------------------------------------------------------------------------------
context('Some tests for imputex')

test_that('Test that return is as expected', {
  expect_is(impute$imputations , 'data.frame') 
  expect_identical(impute$imputations, na.omit(impute$imputations))
  expect_equal(ncol(impute$imputations), impute$number_of_imputations+1) # there are m + 1 columns, as final imputed vector is also attached
  expect_is(impute$fulldata , 'data.frame')
  expect_equal(nrow(impute$fulldata), nrow(rd))
  expect_equal(ncol(impute$fulldata), ncol(rd))
  expect_is(impute$impquantiles , 'data.frame')
  expect_equal(nrow(impute$impquantiles), impute$number_of_imputations)
  expect_identical(impute$impquantiles, na.omit(impute$impquantiles))
  expect_true(all(impute$imputevariance > 0))
})

test_that('Test that call stops correctly with the valid error message', {
  rd$indicator = NULL
  expect_error(imputex(xmu_formula = x1 ~ y + x2, data = rd, indicator = "indicator", censtype = "right"), 'indicator must be a column name in data')
})

cens_values <- rd$x1[rd$indicator == 1] # censored values
imps <- impute$imputations$imputedx     # imputed values

test_that('Valid draws, i.e. they are at least drawn from the valid region',{
  expect_true(all(imps >= cens_values)) # x1 defines the lower bound!
  expect_true(all(imps >= cens_values)) # x1 defines the upper bound!
})
####### Inteval censored fehlen! (folgenedes Muster ungefähr)
#expect_true(all(samplecensored(model ,censtype = 'interval', predict.df, rd, censor = "x1")$draw >= Wdat$cens$lower) &  # interval is not implemented yet
# all(samplecensored(...,censtype = 'interval')$draw <= Wdat$cens$upper))
