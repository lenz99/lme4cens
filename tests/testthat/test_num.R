library(lme4cens)

context("Numerical stuff")

MY_TOL <- 0.0001


test_that("log(x+y)", {
  expect_identical(logxpy(-5, -8), log(exp(-5) + exp(-8)) )
  expect_identical(logxpy(-8, -5), log(exp(-5) + exp(-8)) )
  # vectorized
  expect_identical(logxpy(c(-5, -8, -7), c(-8, -8, -8)), log(exp(c(-5, -8, -7)) + exp(-8)) )
})


test_that("log(x-y)", {
  expect_identical(logxmy(-5, -6),
                   log(exp(-5) - exp(-6)))
  # vectorized
  expect_identical(logxmy(c(-3, -6), c(-8, -8)),
                   log(exp(c(-3, -6)) - exp(-8)))
})


MU <- -1.3
SIGMA2 <- 3.1

test_that("Gauss-Hermite Quadrature", {
  # scaled density of standard normal as integrand
  expect_equal(int_gh(f=1), sqrt(pi), tolerance = 1e-7)
  # expectation of X^2, when X ~ N(MU, SIGMA2)
  expect_equal(stats::integrate(f = function(x) sqrt(2 * pi*SIGMA2)**-1 * exp(-(x-MU)^2/(2*SIGMA2)) * x^2, lower = -Inf, upper = Inf)$value,
               sqrt(pi)**-1 * int_gh(f = function(x) (sqrt(2 * SIGMA2) * x + MU)^2, o = 7), tolerance = 1e-5)
})


test_that("Derivatives of normal density", {

  # standard normal density distribution function
  param_dens0 <- list(x =  seq(-2, 2, len = 5))

  ours_dens0 <- do.call(dnorm, param_dens0)
  attr(ours_dens0, "gradient") <- diag( - param_dens0[["x"]] * ours_dens0)

  expect_equal(ours_dens0,
               numericDeriv(quote(dnorm(x)), theta = "x", rho = list2env(param_dens0)),
               tolerance = MY_TOL)

  # normal density as function of x with mean and sd parameters
  param_dens1 <- list(x = 1,
                 mean = seq(-1, 2, len = 7L),
                 sd = seq(.5, 3, len = 7L))

  ours_dens1 <- do.call(dnorm, args = param_dens1)
  attr(ours_dens1, "gradient") <- as.matrix( - (param_dens1[["x"]] - param_dens1[["mean"]]) / param_dens1[["sd"]]^2 * ours_dens1)

  expect_equal(ours_dens1,
               numericDeriv(quote(dnorm(x, mean, sd)), theta = "x", list2env(param_dens1)),
               tolerance = MY_TOL)


  # normal density as function of σ
  param_dens_sd <- list(x = rnorm(5, sd = 2.3),
                        mean = rnorm(5, sd = .3),
                        sd = 1.1)

  ours_dens_sd <- do.call(dnorm, args = param_dens_sd)
  attr(ours_dens_sd, "gradient") <- as.matrix(
    1 / param_dens_sd[["sd"]] * ours_dens_sd * ( (param_dens_sd[["x"]] - param_dens_sd[["mean"]])^2 / param_dens_sd[["sd"]]^2 - 1)
  )

  expect_equal(ours_dens_sd,
               numericDeriv(quote(dnorm(x, mean, sd)), theta = "sd", list2env(param_dens_sd)),
               tolerance = MY_TOL)


  ## product of normal densities: use param_dens_sd's x for product, results in single value
  ours_dens_sd2 <- prod(ours_dens_sd)
  attr(ours_dens_sd2, "gradient") <- ours_dens_sd2 / param_dens_sd[["sd"]] * ( crossprod(param_dens_sd[["x"]] - param_dens_sd[["mean"]]) / param_dens_sd[["sd"]]^2 - length(ours_dens_sd) )

  expect_equal(ours_dens_sd2,
               numericDeriv(quote( (2*pi*sd^2)^-(length(ours_dens_sd)/2) * exp(-1/(2 * sd^2) * sum((param_dens_sd[["x"]] - param_dens_sd[["mean"]])^2))),
                            theta = "sd", list2env(param_dens_sd)),
               tolerance = MY_TOL)

})



test_that("Derivatives of cumulative normal distribution function", {

  # cum. normal density as function of location parameter (occurs in left-censoring)
  # @param q quantile value
  # @param x mean value from design matrix
  # @param b beta parameter value
  # @param s standard deviation
  pnorm_beta <- function(q, b, x, s) pnorm(q = q, mean = b * x, s = s)
  dnorm_beta <- function(q, b, x, s) dnorm(x = q, mean = b * x, s = s)

  # parameters for derivative of b
  param_loc <- list(q=seq(-1, 1, len = 11),
                 b = 1,
                 x=seq(-2, 2, len = 11),
                 s = sqrt(SIGMA2))

  # solution: derivation w.r.t. location parameter b [as on paper]
  ours_loc <- do.call(pnorm_beta, args = param_loc)
  attr(ours_loc, "gradient") <- as.matrix(
    - do.call(dnorm_beta, args = param_loc) * param_loc[["x"]]
  )

  expect_equal(ours_loc,
               numericDeriv(quote(pnorm_beta(q, b, x, s)), theta = "b", list2env(param_loc)),
               tolerance = MY_TOL)



  # derivatives w.r.t. lb (=beta on log-scale)
  pnorm_lbeta <- function(q, lb, x, s) pnorm(q = q, mean = exp(lb) * x, s = s)
  dnorm_lbeta <- function(q, lb, x, s) dnorm(x = q, mean = exp(lb) * x, s = s)

  param_lloc <- list(q = seq(-1, 1, len = 5),
                 lb = -.1,
                 x = rnorm(5, mean = 0, sd = 2.3),
                 s = sqrt(SIGMA2))

  ours_lloc <- do.call(pnorm_lbeta, args = param_lloc)
  attr(ours_lloc, "gradient") <- as.matrix(
    - do.call(dnorm_lbeta, param_lloc) * param_lloc[["x"]] * exp(param_lloc[["lb"]])
  )

  expect_equal(ours_lloc,
               numericDeriv(quote(pnorm_lbeta(q, lb, x, s)), theta = "lb", list2env(param_lloc)),
               tolerance = MY_TOL)




  # parameters for derivative wrt. s
  param_sd <- list(q = seq(-1, 1, len = 11),
                 b = rnorm(11, sd = 3.2),
                 x = seq(-2, 2, len = 11),
                 s = 2)


  ours_sd <- do.call(pnorm_beta, args = param_sd)
  attr(ours_sd, "gradient") <- as.matrix(
    - do.call(dnorm_beta, param_sd) * (param_sd[["q"]] - param_sd[["b"]] * param_sd[["x"]]) / param_sd[["s"]]
  )

  # # visualization of first parameter set:
  # # curve and its derivative at selected point
  # param2_1 <- sapply(param2, head, n=1)
  # f1 <- function(s) pnorm_beta(a = param2_1["a"], b = param2_1["b"], x = param2_1["x"], s = s)
  # grad_f1 <- attr(ours2, "gradient")[1]
  # curve(f1, from = param2_1["s"]-2, to = param2_1["s"]+2)
  #
  # lines(x = c(param2_1["s"]-1L, param2_1["s"]+1L),
  #       y = c(ours2[1L] - grad_f1, ours2[1L] + grad_f1), col = "green", lty = "dotted")

  expect_equal(ours_sd,
               numericDeriv(quote(pnorm_beta(q, x, b, s)), theta = "s", list2env(param_sd)),
               tolerance = MY_TOL)




  # parameters for derivative wrt. lsigma (=SD on log-scale)
  # this is like the previous case but chained with the exponential
  pnorm_lsigma <- function(q, mu, lsigma) pnorm(q = q, mean = mu, sd = exp(lsigma))
  dnorm_lsigma <- function(q, mu, lsigma) dnorm(x = q, mean = mu, sd = exp(lsigma))

  param_lsd <- list(q = seq(-1, 1, len = 3),
                 mu = rnorm(3, sd = 3.2),
                 lsigma = -.1)

  ours_lsd <- do.call(pnorm_lsigma, args = param_lsd)
  attr(ours_lsd, "gradient") <- as.matrix(
    - do.call(dnorm_lsigma, param_lsd) * (param_lsd[["q"]] - param_lsd[["mu"]])
    )

  expect_equal(ours_lsd,
               numericDeriv(quote(pnorm_lsigma(q, mu, lsigma)), theta = "lsigma", list2env(param_lsd)),
               tolerance = MY_TOL)


})


