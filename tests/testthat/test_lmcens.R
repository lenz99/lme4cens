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

  ## likelihood functions
  negLogLikFun_left <- fm.aff.simple.left$negLogLikFun
  negLogLikGradFun_left <- attr(negLogLikFun_left, "grad")
  negLogLikFun2_left <- function(beta0, beta1, lsigma) as.vector(negLogLikFun_left(c(beta0, beta1, lsigma)))

  ## likelihood functions
  negLogLikFun_right <- fm.aff.simple.right$negLogLikFun
  negLogLikGradFun_right <- attr(negLogLikFun_right, "grad")
  negLogLikFun2_right <- function(beta0, beta1, lsigma) as.vector(negLogLikFun_right(c(beta0, beta1, lsigma)))

  negLogLikFun_lr <- fm.aff.simple.lr$negLogLikFun
  negLogLikGradFun_lr <- attr(negLogLikFun_lr, "grad")
  negLogLikFun2_lr <- function(beta0, beta1, lsigma) as.vector(negLogLikFun_lr(c(beta0, beta1, lsigma)))


  # parameter 1 ------
  param1 <- c(-4, 0.5, 2) # intercept, slope and log(σ)

  ours1_left <- as.vector(negLogLikFun_left(param1))
  attr(ours1_left, "gradient") <- matrix(as.vector(negLogLikGradFun_left(param1)), nrow = 1L)

  ours1_right <- as.vector(negLogLikFun_right(param1))
  attr(ours1_right, "gradient") <- matrix(as.vector(negLogLikGradFun_right(param1)), nrow = 1L)

  ours1_lr <- as.vector(negLogLikFun_lr(param1))
  attr(ours1_lr, "gradient") <- matrix(as.vector(negLogLikGradFun_lr(param1)), nrow = 1L)

  # parameter 2 ------
  param2 <- c(-7, 0.75, 4)

  ours2_left <- as.vector(negLogLikFun_left(param2))
  attr(ours2_left, "gradient") <- matrix(as.vector(negLogLikGradFun_left(param2)), nrow = 1L)

  ours2_right <- as.vector(negLogLikFun_right(param2))
  attr(ours2_right, "gradient") <- matrix(as.vector(negLogLikGradFun_right(param2)), nrow = 1L)

  ours2_lr <- as.vector(negLogLikFun_lr(param2))
  attr(ours2_lr, "gradient") <- matrix(as.vector(negLogLikGradFun_lr(param2)), nrow = 1L)


  setEnv_simple <- function(param){
    evalEnv <- new.env()
    assign("beta0", param[1L], envir = evalEnv)
    assign("beta1", param[2L], envir = evalEnv)
    assign("lsigma", param[3L], envir = evalEnv)

    evalEnv
  }

  expect_equal(ours1_left, numericDeriv(quote(negLogLikFun2_left(beta0, beta1, lsigma)), c("beta0", "beta1", "lsigma"), setEnv_simple(param1)), tolerance = MY_TOL)
  expect_equal(ours2_left, numericDeriv(quote(negLogLikFun2_left(beta0, beta1, lsigma)), c("beta0", "beta1", "lsigma"), setEnv_simple(param2)), tolerance = MY_TOL)

  expect_equal(ours1_right, numericDeriv(quote(negLogLikFun2_right(beta0, beta1, lsigma)), c("beta0", "beta1", "lsigma"), setEnv_simple(param1)), tolerance = MY_TOL)
  expect_equal(ours2_right, numericDeriv(quote(negLogLikFun2_right(beta0, beta1, lsigma)), c("beta0", "beta1", "lsigma"), setEnv_simple(param2)), tolerance = MY_TOL)

  expect_equal(ours1_lr, numericDeriv(quote(negLogLikFun2_lr(beta0, beta1, lsigma)), c("beta0", "beta1", "lsigma"), setEnv_simple(param1)), tolerance = MY_TOL)
  expect_equal(ours2_lr, numericDeriv(quote(negLogLikFun2_lr(beta0, beta1, lsigma)), c("beta0", "beta1", "lsigma"), setEnv_simple(param2)), tolerance = MY_TOL)
})


