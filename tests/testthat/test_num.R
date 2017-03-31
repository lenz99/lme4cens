library(lme4cens)

context("Numerical stuff")

MU <- -1.3
SIGMA2 <- 3.1

test_that("Gauss-Hermite Quadrature", {
  # scaled density of standard normal as integrand
  expect_equal(int_gh(f=1), sqrt(pi), tolerance = 1e-7)
  # expectation of X^2, when X ~ N(MU, SIGMA2)
  expect_equal(stats::integrate(f = function(x) sqrt(2 * pi*SIGMA2)**-1 * exp(-(x-MU)^2/(2*SIGMA2)) * x^2, lower = -Inf, upper = Inf)$value,
               sqrt(pi)**-1 * int_gh(f = function(x) (sqrt(2 * SIGMA2) * x + MU)^2, o = 7), tolerance = 1e-5)
})

