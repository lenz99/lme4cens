library(lme4cens)

context("Basic tests")

test_that("basic arithmetics", {
  expect_identical(1L+1L, 2L)
  expect_equal(sqrt(2L), 1.414214, tolerance = 1e-5)
})
