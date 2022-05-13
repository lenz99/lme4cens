
library('survival')
library('optimx')

context("Simple linear mixed models with censoring: model fitting")

MY_TOL <- 1e-3
MY_TOL_REL <- .0075
MY_SEED <- 12345L
NSIM <- 26L
PROB_CENS <- .1

# random censoring should not bias the fixed effect estimates
test_that("sleepstudy model simulation, random right censoring", code = {

  # start from linear mixed model as ground truth model
  fm <- lme4::lmer(Reaction ~ Days + (1|Subject), data = sleepstudy, REML = FALSE)

  fm.beta <- fixef(fm)
  fm.sigma <- sigma(fm)
  fm.sig01 <- unname(fm.sigma * getME(fm, "theta"))

  # fit a lmercens-model to it
  fmc <- lmercens(Reaction ~ Days + (1|Subject), data = sleepstudy, REML = FALSE,
                  lmerControl(optimizer = "optimx", optCtrl = list(method = "L-BFGS-B")))
  fmc.beta <- fixef(fmc)
  fmc.sigma <- sigma(fmc)
  fmc.sig01 <- sigma(fmc, which = "betw")

  # mkuhn, 2018-02-06: I would like to use bootMer, but I need random censoring in the response prior to refit
  # simulate response variable
  fm.sim0 <- simulate(fm, nsim = NSIM, seed = MY_SEED)

  # add random right censoring, independent of the simulated value
  cens.stat <- matrix(data = 1L-rbinom(n = NROW(sleepstudy) * NSIM, size = 1L, prob = PROB_CENS), ncol = NSIM)


  # fit models on simulated data
  simRes <- parallel::mclapply(X = 1:NSIM, FUN = function(i){
    fmc_sim <- refit(fmc, newresp = Surv(fm.sim0[[i]], event = cens.stat[,i]))
    c(fixef = fixef(fmc_sim), sig01 = sigma(fmc_sim, which = "betw"), sigma = sigma(fmc_sim))
  })

  #  evaluate parameter fits
  simRes_m <- colMeans(matrix(unlist(simRes[!sapply(simRes, inherits, what = "try-error")]), ncol = 4L, byrow = TRUE,
                              dimnames = list(NULL, cols = c("fixef_int", "fixef_days", "sig01", "sigma"))))


  expect_equal(1+1, 2L, tolerance = MY_TOL)
  expect_equivalent(simRes_m["fixef_int"], expected = fmc.beta[1L],  tolerance = abs(fmc.beta[1L]) * MY_TOL_REL)
  expect_equivalent(simRes_m["fixef_days"], expected = fmc.beta[2L],  tolerance = abs(fmc.beta[2L]) * MY_TOL_REL)
  expect_equivalent(simRes_m["sig01"], expected = fmc.sig01,  tolerance = abs(fmc.sig01) * MY_TOL_REL)
  expect_equivalent(simRes_m["sigma"], expected = fmc.sigma,  tolerance = abs(fmc.sigma) * MY_TOL_REL)
})


