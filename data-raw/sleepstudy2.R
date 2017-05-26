# save subset of sleepstudy data

set.seed(12345L)

library(dplyr)


data(sleepstudy, package="lme4")


## add censored status variable (for left- and right-censoring)
# event2: random left- and right censoring
# event3: fixed censoring depending on threshold values for `Reaction`

# coding for status variable:
## 0 right censored
## 1 event observed
## 2 left censored

REACT_L2 <- 300
REACT_R2 <- 310

REACT_L3 <- 212
REACT_R3 <- 350

sleepstudy <- sleepstudy %>%
  mutate_( event2 = ~ if_else(Reaction < REACT_L2 & runif(NROW(.)) < .2,
                              true = 2L, false = if_else(Reaction > REACT_R2 & runif(NROW(.)) < .23, true = 0L, false = 1L)),
           event3 = ~ if_else(Reaction < REACT_L3, true = 2L, false = if_else(Reaction > REACT_R3, true = 0L, false = 1L)),
           Reaction3 = ~ pmin(pmax(Reaction, REACT_L3), REACT_R3)
  )

attr(sleepstudy, "left2") <- REACT_L2
attr(sleepstudy, "right2") <- REACT_R2


attr(sleepstudy, "left3") <- REACT_L3
attr(sleepstudy, "right3") <- REACT_R3


sleepstudy2 <- sleepstudy %>%
  filter_( ~ Subject %in% c(309, 331, 333, 372) ) %>%
  droplevels


save(sleepstudy, file = "data/sleepstudy.rda")
save(sleepstudy2, file = "data/sleepstudy2.rda")

