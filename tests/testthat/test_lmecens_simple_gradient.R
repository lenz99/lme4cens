library(survival)

context("Simple linear mixed models with censoring: gradient of negative log-likelihood function with numeric GH-integration")

MY_TOL <- 1e-4

set.seed(12345L)
# add random censoring on 20% of observations: FALSE indicates censoring
sleepstudy2$event <- sample(c(TRUE, TRUE, TRUE, TRUE, FALSE), size = NROW(sleepstudy2), replace = TRUE)
# for interval censoring (no left censoring as we would need to have -Inf as first Reaction-time)
sleepstudy2$Reaction2 <- sleepstudy2$Reaction + sample(c(0, 0, 0, 0, 0, 0, 0, 3, 5, 8, 8, Inf),
                                                       size = NROW(sleepstudy2), replace = TRUE)

# parameter vector theta: log(σ_b) and log(σ), then intercept, slope,
param1 <- list(lBetwSD = 3.4, lResSD=2.84, beta0=230, beta1=7)
param2 <- list(lBetwSD = 3.84, lResSD=2.91, beta0=260, beta1=6.5)


test_that("no censoring", {


  # build models ------------------------------------------------------------

  # no censoring, only observations
  fm.lmercens_obs <- lmercens(Surv(Reaction) ~ Days + (1|Subject), data = sleepstudy2,
                          REML = FALSE)

  # negative log-likelihood -------------------------------------------------

  negLogLikFun_obs <- fm.lmercens_obs[["negLogLikFun"]]
  negLogLikFun2_obs <- function(lBetwSD, lResSD, beta0, beta1) as.vector(negLogLikFun_obs(c(lBetwSD, lResSD, beta0, beta1)))
  negLogLikGradFun_obs <- attr(negLogLikFun_obs, "grad")


  # build gradient ------

  ours1_obs <- as.vector(negLogLikFun_obs(unlist(param1)))
  attr(ours1_obs, "gradient") <- matrix(as.vector(negLogLikGradFun_obs(unlist(param1))), nrow = 1L)



  ours2_obs <- as.vector(negLogLikFun_obs(unlist(param2)))
  attr(ours2_obs, "gradient") <- matrix(as.vector(negLogLikGradFun_obs(unlist(param2))), nrow = 1L)


  # assess gradient ---------------------------------------------------------

  # at parameter vector  1
  expect_equal(ours1_obs, numericDeriv(quote(negLogLikFun2_obs(lBetwSD, lResSD, beta0, beta1)),
                                       c("lBetwSD", "lResSD", "beta0", "beta1"),
                                       list2env(param1)),
               tolerance = MY_TOL)

  # at parameter vector 2
  expect_equal(ours2_obs, numericDeriv(quote(negLogLikFun2_obs(lBetwSD, lResSD, beta0, beta1)),
                                       c("lBetwSD", "lResSD", "beta0", "beta1"),
                                       list2env(param2)),
               tolerance = MY_TOL)
})






test_that("left censoring", {

  # with left censorings
  fm.lmercens_left <- lmercens(Surv(Reaction, event, type = "left") ~ Days + (1|Subject), data = sleepstudy2,
                              REML = FALSE)


  negLogLikFun_left <- fm.lmercens_left[["negLogLikFun"]]
  negLogLikFun2_left <- function(lBetwSD, lResSD, beta0, beta1) as.vector(negLogLikFun_left(c(lBetwSD, lResSD, beta0, beta1)))
  negLogLikGradFun_left <- attr(negLogLikFun_left, "grad")



  ours1_left <- as.vector(negLogLikFun_left(unlist(param1)))
  attr(ours1_left, "gradient") <- matrix(as.vector(negLogLikGradFun_left(unlist(param1))), nrow = 1L)

  ours2_left <- as.vector(negLogLikFun_left(unlist(param2)))
  attr(ours2_left, "gradient") <- matrix(as.vector(negLogLikGradFun_left(unlist(param2))), nrow = 1L)



  expect_equal(ours1_left, numericDeriv(quote(negLogLikFun2_left(lBetwSD, lResSD, beta0, beta1)),
                                        c("lBetwSD", "lResSD", "beta0", "beta1"),
                                       list2env(param1)),
               tolerance = MY_TOL)


  expect_equal(ours2_left, numericDeriv(quote(negLogLikFun2_left(lBetwSD, lResSD, beta0, beta1)),
                                        c("lBetwSD", "lResSD", "beta0", "beta1"),
                                        list2env(param2)),
               tolerance = MY_TOL)

})



test_that("right censoring", {
  # with right censorings
  fm.lmercens_right <- lmercens(Surv(Reaction, event, type = "right") ~ Days + (1|Subject), data = sleepstudy2,
                                REML = FALSE)

  negLogLikFun_right <- fm.lmercens_right[["negLogLikFun"]]
  negLogLikFun2_right <- function(lBetwSD, lResSD, beta0, beta1) as.vector(negLogLikFun_right(c(lBetwSD, lResSD, beta0, beta1)))
  negLogLikGradFun_right <- attr(negLogLikFun_right, "grad")

  ours1_right <- as.vector(negLogLikFun_right(unlist(param1)))
  attr(ours1_right, "gradient") <- matrix(as.vector(negLogLikGradFun_right(unlist(param1))), nrow = 1L)

  ours2_right <- as.vector(negLogLikFun_right(unlist(param2)))
  attr(ours2_right, "gradient") <- matrix(as.vector(negLogLikGradFun_right(unlist(param2))), nrow = 1L)


  expect_equal(ours1_right, numericDeriv(quote(negLogLikFun2_right(lBetwSD, lResSD, beta0, beta1)),
                                         c("lBetwSD", "lResSD", "beta0", "beta1"),
                                         list2env(param1)),
               tolerance = MY_TOL)

  expect_equal(ours2_right, numericDeriv(quote(negLogLikFun2_right(lBetwSD, lResSD, beta0, beta1)),
                                         c("lBetwSD", "lResSD", "beta0", "beta1"),
                                         list2env(param2)),
               tolerance = MY_TOL)

})


test_that("interval censoring", {
  # with interval and right censorings
  fm.lmercens_int <- lmercens(Surv(Reaction, time2 = Reaction2, type = "interval2") ~ Days + (1|Subject), data = sleepstudy2,
                              REML = FALSE)

  negLogLikFun_int <- fm.lmercens_int[["negLogLikFun"]]
  negLogLikFun2_int <- function(lBetwSD, lResSD, beta0, beta1) as.vector(negLogLikFun_int(c(lBetwSD, lResSD, beta0, beta1)))
  negLogLikGradFun_int <- attr(negLogLikFun_int, "grad")

  ours1_int <- as.vector(negLogLikFun_int(unlist(param1)))
  attr(ours1_int, "gradient") <- matrix(as.vector(negLogLikGradFun_int(unlist(param1))), nrow = 1L)

  ours2_int <- as.vector(negLogLikFun_int(unlist(param2)))
  attr(ours2_int, "gradient") <- matrix(as.vector(negLogLikGradFun_int(unlist(param2))), nrow = 1L)

  expect_equal(ours1_int, numericDeriv(quote(negLogLikFun2_int(lBetwSD, lResSD, beta0, beta1)),
                                       c("lBetwSD", "lResSD", "beta0", "beta1"),
                                       list2env(param1)),
               tolerance = MY_TOL)


  expect_equal(ours2_int, numericDeriv(quote(negLogLikFun2_int(lBetwSD, lResSD, beta0, beta1)),
                                       c("lBetwSD", "lResSD", "beta0", "beta1"),
                                       list2env(param2)),
               tolerance = MY_TOL)
})


