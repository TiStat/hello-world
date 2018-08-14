library(testthat)

context('evaluate family functions')
test_that('Mismatchning arguments are stoped',{
  expect_error(family_fun(..., func = 'r', ..., p = c(...)), "One of x,q,p,n,... arguments doesn't match with the distributional function (e.g. dNO, pNO, qNO, rNO). See the family's documentation for admissable arguments.")
  expect_is(family_fun(..., func = 'r', n = 10), 'vector')
  expect_is(family_fun(..., func = 'p', q = c(...)), 'vector')
  expect_equal(length(family_fun(..., func = 'r', n = 10)), 10)
})

test_that('return is vector',{
  expect_is(family_fun(...), 'vector')
  expect_equal(nrow(predictdata),length(family_fun(...))) # predictdata muss hier also in der enviroment sein um getestet werden zu können
  expect_identical(family_fun(...), na.omit(family_fun(...)))
})

context('evaluate inverse sampling in samplecensored')
test_that('valid draws. They are at least in drawn from the valid region',{
  expect_success(all(samplecensored(...,censtype = 'right')$draw >= Wdat$cens$x1)) # x1 defines the lower bound!
  expect_success(all(samplecensored(...,censtype = 'left')$draw <= Wdat$cens$x1))
  expect_success(all(samplecensored(...,censtype = 'interval')$draw >= Wdat$cens$lower) &  # interval is not implemented yet
                   all(samplecensored(...,censtype = 'interval')$draw <= Wdat$cens$upper))
})
test_that('Valid quantiles; the distribution\'s quantiles fullfill the constraint of the censored variable', {
  expect_is(samplecensored(...)$quantiles, 'data.frame')
  expect_equal(nrow(samplecensored(...)$quantiles), nrow(predictdata)) # predictdata muss in enviroment sein
  expect_success(all(samplecensored(..., censtype = 'right')$quantiles > Wdat$cens$x1)) # strict, as 0.05 is lowest # comparison is unambigous! vector is recycled
  expect_success(all(samplecensored(..., censtype = 'left')$quantiles < Wdat$cens$x1)) # strict, as 0.95 is highest
  expect_success(all(samplecensored(..., censtype = 'interval')$quantiles > Wdat$cens$lower)) # strict, as 0.05 is lowest
  expect_success(all(samplecensored(..., censtype = 'interval')$quantiles > Wdat$cens$upper)) # strict, as 0.05 is lowest
})