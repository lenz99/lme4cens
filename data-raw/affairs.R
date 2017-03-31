library(dplyr)


data("Affairs", package = "AER")

Affairs %>%
  dplyr::mutate_(event =  ~dplyr::if_else(affairs > 0, true = 1L, false = 0L), # 1=event, 0= left censored
                 # event coding for interval-censoring, type 'interval' (not 'interval2')
                 event2 = ~dplyr::if_else(affairs <= 0, true = 2L, ## 2 = left-censored
                                          false = dplyr::if_else(affairs > 10L, true = 0L, ## 0 = right-censored,
                                                                 false = 1L))) ->
  Affairs


save(Affairs, file = "data/Affairs.rda")
