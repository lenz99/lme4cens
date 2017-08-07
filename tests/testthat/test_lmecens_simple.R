
context("Simple linear mixed models with censoring")

MY_TOL <- 1e-4

test_that("Gradient of lmecens ML function with numeric GH-integration for simple linear mixed model", {

  set.seed(12345L)

  # add random censoring
  sleepstudy2$event <- sample(c(TRUE, TRUE, TRUE, TRUE, FALSE), size = NROW(sleepstudy2), replace = TRUE)
  sleepstudy2$Reaction2 <- sleepstudy2$Reaction + sample(c(0, 0, 0, 0, 0, 0, 0, 3, 5, 8, 8, Inf), size = NROW(sleepstudy2), replace = TRUE)

  paramS <- c(260, 5, 3, 2)


  # build models ------------------------------------------------------------

  # no censoring, only observations
  fm.lmercens.obs <- lmercens(Surv(Reaction) ~ Days + (1|Subject), data = sleepstudy2,
                          start = paramS, REML = FALSE)


  # with left censorings
  fm.lmercens.left <- lmercens(Surv(Reaction, event, type = "left") ~ Days + (1|Subject), data = sleepstudy2,
                              start = paramS, REML = FALSE)

  # with right censorings
  fm.lmercens.right <- lmercens(Surv(Reaction, event, type = "right") ~ Days + (1|Subject), data = sleepstudy2,
                               start = paramS, REML = FALSE)


  # with interval (and right) censorings
  fm.lmercens.int <- lmercens(Surv(Reaction, time2 = Reaction2, type = "interval2") ~ Days + (1|Subject), data = sleepstudy2,
                                start = paramS, REML = FALSE)


  # negative log-likelihood -------------------------------------------------

  negLogLikFun_obs <- fm.lmercens.obs[["negLogLikFun"]]
  negLogLikGradFun_obs <- attr(negLogLikFun_obs, "grad")
  negLogLikFun2_obs <- function(beta0, beta1, lBetwSD, lResSD) as.vector(negLogLikFun_obs(c(beta0, beta1, lBetwSD, lResSD)))


  negLogLikFun_left <- fm.lmercens.left[["negLogLikFun"]]
  negLogLikGradFun_left <- attr(negLogLikFun_left, "grad")
  negLogLikFun2_left <- function(beta0, beta1, lBetwSD, lResSD) as.vector(negLogLikFun_left(c(beta0, beta1, lBetwSD, lResSD)))


  negLogLikFun_right <- fm.lmercens.right[["negLogLikFun"]]
  negLogLikGradFun_right <- attr(negLogLikFun_right, "grad")
  negLogLikFun2_right <- function(beta0, beta1, lBetwSD, lResSD) as.vector(negLogLikFun_right(c(beta0, beta1, lBetwSD, lResSD)))

  negLogLikFun_int <- fm.lmercens.int[["negLogLikFun"]]
  negLogLikGradFun_int <- attr(negLogLikFun_int, "grad")
  negLogLikFun2_int <- function(beta0, beta1, lBetwSD, lResSD) as.vector(negLogLikFun_int(c(beta0, beta1, lBetwSD, lResSD)))



  # build gradient ------

  param1 <- list(beta0=230, beta1=7, lBetwSD = 3.4, lResSD=2.84) # intercept, slope, log(σ_b) and log(σ)

  ours1_obs <- as.vector(negLogLikFun_obs(unlist(param1)))
  attr(ours1_obs, "gradient") <- matrix(as.vector(negLogLikGradFun_obs(unlist(param1))), nrow = 1L)

  ours1_left <- as.vector(negLogLikFun_left(unlist(param1)))
  attr(ours1_left, "gradient") <- matrix(as.vector(negLogLikGradFun_left(unlist(param1))), nrow = 1L)

  ours1_right <- as.vector(negLogLikFun_right(unlist(param1)))
  attr(ours1_right, "gradient") <- matrix(as.vector(negLogLikGradFun_right(unlist(param1))), nrow = 1L)

  ours1_int <- as.vector(negLogLikFun_int(unlist(param1)))
  attr(ours1_int, "gradient") <- matrix(as.vector(negLogLikGradFun_int(unlist(param1))), nrow = 1L)


  # assess gradient ---------------------------------------------------------

  expect_equal(ours1_obs, numericDeriv(quote(negLogLikFun2_obs(beta0, beta1, lBetwSD, lResSD)),
                                        c("beta0", "beta1", "lBetwSD", "lResSD"),
                                       list2env(param1)),
               tolerance = MY_TOL)

  expect_equal(ours1_left, numericDeriv(quote(negLogLikFun2_left(beta0, beta1, lBetwSD, lResSD)),
                                       c("beta0", "beta1", "lBetwSD", "lResSD"),
                                       list2env(param1)),
               tolerance = MY_TOL)


  expect_equal(ours1_right, numericDeriv(quote(negLogLikFun2_right(beta0, beta1, lBetwSD, lResSD)),
                                        c("beta0", "beta1", "lBetwSD", "lResSD"),
                                        list2env(param1)),
               tolerance = MY_TOL)

  expect_equal(ours1_int, numericDeriv(quote(negLogLikFun2_int(beta0, beta1, lBetwSD, lResSD)),
                                         c("beta0", "beta1", "lBetwSD", "lResSD"),
                                         list2env(param1)),
               tolerance = MY_TOL)

  expect_equal(1+1, 2L)
})
