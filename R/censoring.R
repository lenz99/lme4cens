## Handle Censoring


#' Prepare the survival response.
#' It uses the interval-coding that supports left-, right- and interval-censoring.
#' @param y0 response vector of class [survival::Surv].
#' @return response vector, normalized to interval-coding.
prepSurvResp <- function(y0) {
  stopifnot( inherits(y0, what = "Surv") )

  survType <- attr(y0, which = "type")

  if (survType == 'interval') y0 else {
    yMat <- as.matrix(y0)
    yTime <- yMat[, "time"]
    yStat <- yMat[, "status"]
    censLevel <- if (max(yStat) == 2) 1 else 0
    isCens <- (yStat == censLevel)

    switch(survType,
           right={Surv(time  = yTime,
                       time2 = ifelse(isCens, yes = NA_real_, no = yTime),
                       type = "interval2")},
           left={ Surv(time  = ifelse(isCens, yes = NA_real_, no = yTime),
                       time2 = yTime,
                       type = "interval2")},
           stop("this type of censoring is not supported!")
    )
  }
}



#' Get an quick approximation for a flat response vector from the Surv-response.
#' It is useful to get sensible start-values using ordinary regression.
flattenResponse <- function(yTime1, yTime2, yStat) {
  flat_y <- yTime1
  # interval cens: use mid-point
  flat_y[yStat == 3L] <- (flat_y[yStat == 3L] + yTime2[yStat == 3L]) / 2L

  flat_y_sd <- max(1.05 * sd(flat_y), 1.05 * mad(flat_y), .001,
                   na.rm = TRUE)

  # right cens
  flat_y[yStat == 0L] <- flat_y[yStat == 0L] + abs(rnorm(n = sum(yStat == 0L), mean = 0L, sd = 2*flat_y_sd))
  # left cens
  flat_y[yStat == 2L] <- flat_y[yStat == 2L] - abs(rnorm(n = sum(yStat == 2L), mean = 0L, sd = 2*flat_y_sd))

  flat_y
}

