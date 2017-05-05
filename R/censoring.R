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
