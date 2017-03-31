# mkuhn, 2017-02-24
# internal coding of Surv-objects

library(survival)

obs <- data.frame(t1=c(NA_real_, 3, 4, 7),
           t2=c(5, 3, NA_real_, 9))

su <- Surv(time = obs$t1, time2= obs$t2, type = "interval2")

as.matrix(Surv(time = ifelse(is.na(obs$t1), obs$t2, obs$t1), !is.na(obs$t1), type = "left"))
