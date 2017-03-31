## ----setup, include=FALSE, message=FALSE---------------------------------
library(dplyr)
devtools::dev_mode(on = TRUE)
## knitr package options
knitr::opts_knit$set(verbose = FALSE)
options(digits = 3L)

library(survival)
library(lme4cens)
library(microbenchmark)

data("Affairs", package = "lme4cens")

VERBOSE <- 0L

## ----expLeftCens---------------------------------------------------------
set.seed(123L)
N <- 230L
λ <- 1.8

ex_exp_lCens <- data_frame_(list(
  X = ~ rweibull(n = N, shape = 1L, scale = 1/λ),
  Cl = ~ rweibull(n = N, shape = 1.2, scale = 1/λ-.3), # shorter mean time than for response X
  T = ~ pmax(X, Cl),
  status = ~ as.numeric(X >= Cl)
))


log1mexp <- function(x){
  stopifnot( x < 0L, is.numeric(x) )
  
  log1pBranch <- x < -log(2)
  ##y <- numeric(length(x))
  x[log1pBranch] <- log1p(-exp(x[log1pBranch]))
  x[!log1pBranch] <- log(-expm1(x[!log1pBranch]))
  
  x
}

# factory method for negative log-likelihood function
negLogLik <- function(data){
  stopifnot( is.data.frame(data), all(c("T", "status") %in% names(data)) )
  
  nbrEvents <- sum(data$status)
  eventIdx <- data$status == 1L
  
  function(l)
    - nbrEvents * log(l) + l * sum(data$T[eventIdx]) - sum(log1mexp(-l * data$T[!eventIdx]))
}


optimize(negLogLik(data = ex_exp_lCens),
         interval = c(0, 10))

## ----expLeftCens_survreg-------------------------------------------------
survreg(Surv(T, status, type = "left") ~ 1, dist = "exp",
        data = ex_exp_lCens)

## ----ex1_lmcens_affairs--------------------------------------------------
fm <- lmcens(Surv(affairs, event, type = "left") ~ age + yearsmarried + religiousness + occupation + rating, data = Affairs)

summary(fm)

## ----ex1_lmcens_affairs_w------------------------------------------------

wPos <- with(Affairs, which(affairs > 1 & religiousness > 2))
w <- rep(1L, 601)
w[wPos] <- 2.5


fmw <- lmcens(Surv(affairs, event, type = "left") ~ age + yearsmarried + religiousness + occupation + rating,
              weights = w, data = Affairs)
summary(fmw)

