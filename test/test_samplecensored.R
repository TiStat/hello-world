library(testthat)

context('evaluate family functions')
test_that('Mismatchning arguments are stoped',{
  expect_error(family_fun(object = , func = 'r', predictdata = , fitdata = , p = c()), "One of x,q,p,n,... arguments doesn't match with the distributional function (e.g. dNO, pNO, qNO, rNO). See the family's documentation for admissable arguments.")
})

test_that('return is vector',{
  expect_is(family_fun(...), 'vector')
  expect_equal(nrow(predictdata),length(family_fun(...))) # predictdata muss hier also in der enviroment sein um getestet werden zu können
  expect_identical(family_fun(...), na.omit(family_fun(...)))
})


context('evaluate inverse sampling in samplecensored')
test_that('valid draws',{
  expect_success(all(samplecensored(...,censtype = 'right')$draw >= Wdat$cens$x1))
  expect_success(all(samplecensored(...,censtype = 'left')$draw <= Wdat$cens$x1))
  expect_success(all(samplecensored(...,censtype = 'interval')$draw >= Wdat$cens$lower) & 
                   all(samplecensored(...,censtype = 'interval')$draw <= Wdat$cens$upper))
})
test_that('quantiles', {
  expect_is(samplecensored(...)$quantiles, 'data.frame')
  expect_equal(nrow(samplecensored(...)$quantiles), nrow(predictdata))
})