library(testthat)

context('evaluate family functions')
test_that('Mismatchning arguments are stoped',{
  expect_error(family_fun(..., func = 'r', ..., p = c(...)), "One of x,q,p,n,... arguments doesn't match with the distributional function (e.g. dNO, pNO, qNO, rNO). See the family's documentation for admissable arguments.")
})

test_that('output as expected', {
  expect_is(family_fun(..., func = 'r', n = 10), 'vector')
  expect_is(family_fun(..., func = 'p', q = c(...)), 'vector')
  expect_equal(length(family_fun(..., func = 'r', n = 10)), 10)

  # next is only valid, if func != 'r'
  expect_equal(nrow(predictdata),length(family_fun(...))) # predictdata muss hier also in der enviroment sein um getestet werden zu können
  expect_identical(family_fun(...), na.omit(family_fun(...)))
})

context('evaluate inverse sampling in samplecensored')
test_that('valid draws. They are at least in drawn from the valid region',{
  expect_true(all(samplecensored(...,censtype = 'right')$draw >= Wdat$cens$x1)) # x1 defines the lower bound!
  expect_true(all(samplecensored(...,censtype = 'left')$draw <= Wdat$cens$x1))
  expect_true(all(samplecensored(...,censtype = 'interval')$draw >= Wdat$cens$lower) &  # interval is not implemented yet
                   all(samplecensored(...,censtype = 'interval')$draw <= Wdat$cens$upper))
})
test_that('Valid quantiles; the distribution\'s quantiles fullfill the constraint of the censored variable', {
  expect_is(samplecensored(...)$quantiles, 'data.frame')
  expect_equal(nrow(samplecensored(...)$quantiles), nrow(predictdata)) # predictdata muss in enviroment sein
  expect_true(all(samplecensored(..., censtype = 'right')$quantiles > Wdat$cens$x1)) # strict, as 0.05 is lowest # comparison is unambigous! vector is recycled
  expect_true(all(samplecensored(..., censtype = 'left')$quantiles < Wdat$cens$x1)) # strict, as 0.95 is highest
  expect_true(all(samplecensored(..., censtype = 'interval')$quantiles > Wdat$cens$lower)) # strict, as 0.05 is lowest
  expect_true(all(samplecensored(..., censtype = 'interval')$quantiles > Wdat$cens$upper)) # strict, as 0.05 is lowest
})

data = simulateData(...)
impute = imputex(...)

context('imputex')
test_that('return is as expected', {
  expect_is(impute$imputemat , 'data.frame')
  expect_identical(impute$imputemat, na.omit(impute$imputemat))
  expect_equal(ncol(impute$imputemat), impute$number_of_imputations+1) # there are m + 1 columns, as final imputed vector is also attached
  expect_is(impute$fulldata , 'data.frame')
  expect_equal(nrow(impute$fulldata), nrow(data))
  expect_equal(ncol(impute$fulldata), ncol(data))
  expect_is(impute$impquantiles , 'data.frame')
  expect_equal(nrow(impute$impquantiles), impute$number_of_imputations)
  expect_identical(impute$impquantiles, na.omit(impute$impquantiles))
  expect_true(impute$imputevariance > 0)
})

data = simulateData(...)
data$indicator = NULL
test_that('stops correctly', {
  expect_error(imputex(..., data), 'indicator must be a column name in data')
})

