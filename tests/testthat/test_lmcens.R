library(lme4cens)

context("Linear model with censoring [lmcens] stuff")

MY_TOL <- 1e-3

test_that("Gradient of lmcens ML function: simple regression with different censorings", {


  # fit a simple regression with left-censored response
  fm.aff.simple.left <- lmcens(survival::Surv(affairs, event, type = "left") ~ yearsmarried, data = Affairs)
  # fit a simple regression with right-censored response: actually saying  0 means "at least 0"
  fm.aff.simple.right <- lmcens(survival::Surv(affairs, event, type = "right") ~ yearsmarried, data = Affairs)
  # model with left and right censoring
  fm.aff.simple.lr <- lmcens(survival::Surv(time = affairs, time2 = NULL, event = event2, type = "interval") ~ yearsmarried, data = Affairs)
  # model with simulated interval censoring
  set.seed(12345L)
  Affairs$affairs1 <- Affairs$affairs
  Affairs$affairs1[Affairs$affairs1 == 0L & rbinom(n = NROW(Affairs), size = 1L, prob = .2) == 1L] <- -Inf
  Affairs$affairs2 <- Affairs$affairs + is.finite(Affairs$affairs1) * rbinom(n = NROW(Affairs), size = 3L, prob = .15)
  Affairs$affairs2[Affairs$affairs2 > 0L  & rbinom(n = NROW(Affairs), size = 1L, prob = .2) == 1L] <- +Inf
  fm.aff.simple.int <- lmcens(survival::Surv(time = affairs1, time2 = affairs2, type = "interval2") ~ yearsmarried, data = Affairs)


  ## likelihood functions for left-cens
  negLogLikFun_left <- fm.aff.simple.left$negLogLikFun
  negLogLikGradFun_left <- attr(negLogLikFun_left, "grad")
  negLogLikFun2_left <- function(beta0, beta1, lsigma) as.vector(negLogLikFun_left(c(beta0, beta1, lsigma)))

  ## likelihood functions for right-cens
  negLogLikFun_right <- fm.aff.simple.right$negLogLikFun
  negLogLikGradFun_right <- attr(negLogLikFun_right, "grad")
  negLogLikFun2_right <- function(beta0, beta1, lsigma) as.vector(negLogLikFun_right(c(beta0, beta1, lsigma)))

  ## likelihood functions for left-right-cens
  negLogLikFun_lr <- fm.aff.simple.lr$negLogLikFun
  negLogLikGradFun_lr <- attr(negLogLikFun_lr, "grad")
  negLogLikFun2_lr <- function(beta0, beta1, lsigma) as.vector(negLogLikFun_lr(c(beta0, beta1, lsigma)))

  ## likelihood for interval-cens
  negLogLikFun_int <- fm.aff.simple.int$negLogLikFun
  negLogLikGradFun_int <- attr(negLogLikFun_int, "grad")
  negLogLikFun2_int <- function(beta0, beta1, lsigma) as.vector(negLogLikFun_int(c(beta0, beta1, lsigma)))



  # parameter 1 ------

  param1 <- list(beta0=-4, beta1=0.5, lsigma=2) # intercept, slope and log(σ)

  ours1_left <- as.vector(negLogLikFun_left(unlist(param1)))
  attr(ours1_left, "gradient") <- matrix(as.vector(negLogLikGradFun_left(unlist(param1))), nrow = 1L)

  ours1_right <- as.vector(negLogLikFun_right(unlist(param1)))
  attr(ours1_right, "gradient") <- matrix(as.vector(negLogLikGradFun_right(unlist(param1))), nrow = 1L)

  ours1_lr <- as.vector(negLogLikFun_lr(unlist(param1)))
  attr(ours1_lr, "gradient") <- matrix(as.vector(negLogLikGradFun_lr(unlist(param1))), nrow = 1L)

  ours1_int <- as.vector(negLogLikFun_int(unlist(param1)))
  attr(ours1_int, "gradient") <- matrix(as.vector(negLogLikGradFun_int(unlist(param1))), nrow = 1L)



  # parameter 2 ------

  param2 <- list(beta0=-7, beta1= 0.75, lsigma = 4)

  ours2_left <- as.vector(negLogLikFun_left(unlist(param2)))
  attr(ours2_left, "gradient") <- matrix(as.vector(negLogLikGradFun_left(unlist(param2))), nrow = 1L)

  ours2_right <- as.vector(negLogLikFun_right(unlist(param2)))
  attr(ours2_right, "gradient") <- matrix(as.vector(negLogLikGradFun_right(unlist(param2))), nrow = 1L)

  ours2_lr <- as.vector(negLogLikFun_lr(unlist(param2)))
  attr(ours2_lr, "gradient") <- matrix(as.vector(negLogLikGradFun_lr(unlist(param2))), nrow = 1L)

  ours2_int <- as.vector(negLogLikFun_int(unlist(param2)))
  attr(ours2_int, "gradient") <- matrix(as.vector(negLogLikGradFun_int(unlist(param2))), nrow = 1L)



  expect_equal(ours1_left, numericDeriv(quote(negLogLikFun2_left(beta0, beta1, lsigma)),
                                        c("beta0", "beta1", "lsigma"), list2env(param1)),
               tolerance = MY_TOL)
  expect_equal(ours2_left, numericDeriv(quote(negLogLikFun2_left(beta0, beta1, lsigma)),
                                        c("beta0", "beta1", "lsigma"), list2env(param2)),
               tolerance = MY_TOL)

  expect_equal(ours1_right, numericDeriv(quote(negLogLikFun2_right(beta0, beta1, lsigma)),
                                         c("beta0", "beta1", "lsigma"), list2env(param1)),
               tolerance = MY_TOL)
  expect_equal(ours2_right, numericDeriv(quote(negLogLikFun2_right(beta0, beta1, lsigma)),
                                         c("beta0", "beta1", "lsigma"), list2env(param2)),
               tolerance = MY_TOL)

  expect_equal(ours1_lr, numericDeriv(quote(negLogLikFun2_lr(beta0, beta1, lsigma)),
                                      c("beta0", "beta1", "lsigma"), list2env(param1)),
               tolerance = MY_TOL)
  expect_equal(ours2_lr, numericDeriv(quote(negLogLikFun2_lr(beta0, beta1, lsigma)),
                                      c("beta0", "beta1", "lsigma"), list2env(param2)),
               tolerance = MY_TOL)

  expect_equal(ours1_int, numericDeriv(quote(negLogLikFun2_int(beta0, beta1, lsigma)),
                                       c("beta0", "beta1", "lsigma"), list2env(param1)),
               tolerance = MY_TOL)
  expect_equal(ours2_int, numericDeriv(quote(negLogLikFun2_int(beta0, beta1, lsigma)),
                                       c("beta0", "beta1", "lsigma"), list2env(param2)),
               tolerance = MY_TOL)

})


