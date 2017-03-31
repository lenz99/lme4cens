# save subset of sleepstudy data

set.seed(12345L)

library(dplyr)


data(sleepstudy, package="lme4")


## add censored response (for left- and right-censoring)
## 0 right censored
## 1 event observed
## 2 left censored

# sleepstudy <- sleepstudy %>%
#   mutate_( event2 = ~ if_else(Reaction < 300 & runif(NROW(.)) < .2,
#                               true = 2L, false = if_else(Reaction > 310 & runif(NROW(.)) < .23, true = 0L, false = 1L))
#   )
sleepstudy <- sleepstudy %>%
  mutate_( event2 = ~ if_else(Reaction < 300 & runif(NROW(.)) < .2,
                              true = 2L, false = if_else(Reaction > 310 & runif(NROW(.)) < .23, true = 0L, false = 1L)),
           event3 = ~ if_else(Reaction < 212, true = 2L, false = if_else(Reaction > 350, true = 0L, false = 1L))
  )



sleepstudy2 <- sleepstudy %>%
  filter_( ~ Subject %in% c(309, 331, 333, 372) ) %>%
  droplevels

save(sleepstudy, file = "data/sleepstudy.rda")
save(sleepstudy2, file = "data/sleepstudy2.rda")

