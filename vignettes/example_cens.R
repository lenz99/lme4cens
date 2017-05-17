## ----setup, include=FALSE, message=FALSE---------------------------------
library(dplyr)

devtools::dev_mode(on = TRUE)
## knitr package options
knitr::opts_knit$set(verbose = FALSE)
options(digits = 3L)

library(survival)

library(lme4cens)
library(microbenchmark)

library(censReg)
library(crch)

data("Affairs", package = "lme4cens")

VERBOSE <- 0L

## ----Surv_obj_mat, echo = FALSE------------------------------------------
survObj <- Surv(time = c(3, -Inf, 3, 3), time2 = c(3, 3, Inf, 5), type = "interval2")
survObj
as.matrix(survObj)

## ----expLeftCens_data----------------------------------------------------
set.seed(123L)
N <- 230L
λ <- 1.8

ex_exp_lCens <- data_frame_(list(
  X = ~ rweibull(n = N, shape = 1L, scale = 1/λ),
  Cl = ~ rweibull(n = N, shape = 1.2, scale = 1/λ-.3), # shorter mean time than for response X
  T = ~ pmax(X, Cl),
  status = ~ as.numeric(X >= Cl)
))

## ----expLeftCens_logLik--------------------------------------------------

# factory method for negative log-likelihood function
# function log1mexp efficiently calculates log(1-exp(x))
negLogLik <- function(data){
  stopifnot( is.data.frame(data), all(c("T", "status") %in% names(data)) )
  
  nbrEvents <- sum(data$status)
  eventIdx <- data$status == 1L
  
  function(λ)
    - nbrEvents * log(λ) + λ * sum(data$T[eventIdx]) - sum(lme4cens::log1mexp(-λ * data$T[!eventIdx]))
}


optimize(negLogLik(data = ex_exp_lCens),
         interval = c(0, 10))

## ----expLeftCens_survreg-------------------------------------------------
fm.exp <- survreg(Surv(T, status, type = "left") ~ 1, dist = "exp",
        data = ex_exp_lCens)
fm.exp

## ----affairs_event_boxplot-----------------------------------------------
boxplot(affairs ~ event, data = Affairs, xlab = "Variable 'event'", ylab = "Number of affairs")

## ----affairs_simple_lmcens-----------------------------------------------
fm.aff.simple <- lmcens(Surv(affairs, event, type = "left") ~ yearsmarried, data = Affairs)
summary(fm.aff.simple)

fm.aff.simple.bfgs <- lmcens(Surv(affairs, event, type = "left") ~ yearsmarried,
                             method = "BFGS",
                             data = Affairs)
summary(fm.aff.simple.bfgs)

## ----affairs_simple_censReg, echo = FALSE--------------------------------
fm.simple.aff.censreg <- censReg(affairs ~ yearsmarried, left = 0L, data = Affairs)

summary(fm.simple.aff.censreg)

## ----affairs_multiple_lmcens---------------------------------------------
fm.aff <- lmcens(Surv(affairs, event, type = "left") ~ age + yearsmarried + religiousness + occupation + rating,
                 data = Affairs)

summary(fm.aff)

## ----affairs_multiple_censReg--------------------------------------------
fm.aff.bfgs <- lmcens(Surv(affairs, event, type = "left") ~ age + yearsmarried + religiousness + occupation + rating,
                      method = "BFGS",
                      data = Affairs)

summary(fm.aff.bfgs)

## ----affairs_multiple_censReg__NelderMead--------------------------------
fm.aff.censreg <- censReg(affairs ~ age + yearsmarried + religiousness + occupation + rating, 
                          left = 0L,
                          data = Affairs)

summary(fm.aff.censreg)

## ----ex1_compare_logLiks-------------------------------------------------
llContribs <- attr(fm.aff$logLik.contribs, "loglik.contribs")
llContribs.censreg <- update(fm.aff.censreg, start = coef(fm.aff), logLikOnly = TRUE)

MASS::eqscplot(x = sort(llContribs), y = sort(llContribs.censreg)); abline(0,1, col = "grey", lty = 2)

## ----ex1_lmcens_affairs_censRegStart-------------------------------------
fm.aff2 <- lmcens(Surv(affairs, event, type = "left") ~ age + yearsmarried + religiousness + occupation + rating,
                  start = round(coef(fm.aff.censreg),3),
                  data = Affairs)

summary(fm.aff2)

## ----ex1_compare_logLiks2------------------------------------------------
llContribs2 <- attr(fm.aff2$logLik.contribs, "loglik.contribs")
llContribs2.censreg <- update(fm.aff.censreg, start = coef(fm.aff2), logLikOnly = TRUE)

MASS::eqscplot(x = sort(llContribs2), y = sort(llContribs2.censreg)); abline(0,1, col = "grey", lty = 2)

## ----ex1_crch_affairs----------------------------------------------------

fm.crch.aff <- crch(affairs ~ age + yearsmarried + religiousness + occupation + rating,
                    data = Affairs, left = 0, dist = "gaussian")
summary(fm.crch.aff)

## ----affairs_lmcens_w----------------------------------------------------
wPos <- with(Affairs, which(affairs > 1 & religiousness > 2))
w <- rep(1L, NROW(Affairs))
w[wPos] <- 2.5


fm.aff.w.nm <- lmcens(Surv(affairs, event, type = "left") ~ age + yearsmarried + religiousness + occupation + rating,
              weights = w, data = Affairs)
summary(fm.aff.w.nm)

## ----affairs_lmcens_w_bfgs-----------------------------------------------


fm.aff.w.bfgs <- lmcens(Surv(affairs, event, type = "left") ~ age + yearsmarried + religiousness + occupation + rating,
                        method = "BFGS",
                        weights = w, data = Affairs)
summary(fm.aff.w.bfgs)

## ----ex2_sleepstudy2_lme4cens--------------------------------------------
REACT_L <- 212
REACT_R <- 350

sleepstudy2 %>% 
  lFormula(Surv(pmin(pmax(Reaction, REACT_L), REACT_R), time2 = Reaction, event = event3,
                type = "interval") ~ Days + (1|Subject), data = ., REML=FALSE) -> 
  lForm



lForm %>% 
  append(list(verbose=VERBOSE, quadrature = "stats")) %>%
  do.call(mkLmerCensDevfun_rInt_R, .) ->
  
  myDevFun_f

# Gauss-Hermite quadrature
lForm %>% 
  append(list(verbose=VERBOSE, quadrature = "gh")) %>%
  do.call(mkLmerCensDevfun_rInt_R, .) ->
  
  myDevFun_gh_f


## optimize function
paramS <- c(260, 5, 3, 2)

paramEst <- optim(par = paramS, fn = myDevFun_gh_f)
paramEst

## ----ex2_sleepstudy2_censReg, message = FALSE----------------------------

sleepstudy2 %>% 
  dplyr::mutate_(Reaction = ~ pmin(pmax(Reaction, REACT_L), REACT_R))  %>% 
  plm::pdata.frame(index = c("Subject", "Days")) ->
  
  sleepstudy2.cr

fm.censReg <- censReg::censReg(Reaction ~ as.numeric(Days), left = REACT_L, right = REACT_R,
                               start = paramS, nGHQ = 8L, data = sleepstudy2.cr)

summary(fm.censReg)

## ----sleepstudy2_desc, echo = FALSE, results='asis'----------------------
sleepstudy2 %>%
  dplyr::group_by_(~ Days) %>%
  dplyr::summarise_(mReaction = ~ mean(Reaction)) %>% 
  knitr::kable(digits = 0L)

## ----ex2_sleepstudy2_lme4censLL------------------------------------------
optim(par = coef(fm.censReg), fn = myDevFun_gh_f)

## ----ex2_sleepstudy2_censRegLL, echo = FALSE-----------------------------
ll.censReg <- censReg::censReg(Reaction ~ as.numeric(Days), left = REACT_L, right = REACT_R, start = paramEst$par, logLikOnly = TRUE, nGHQ = 8L,
                               data = sleepstudy2.cr)

nll.lme4cens <- myDevFun_gh_f(param = paramEst$par)
nll.lme4cens.numIntStats <- myDevFun_f(param = paramEst$par)

