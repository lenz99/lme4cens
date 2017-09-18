#' Objective function for LMER-models implemented in R.
#'
#' This implements the objective function (deviance or restricted deviance).
#' This is a test case to compare with the official lme4-implementation.
#' ZZZ offset?
#'
#' @param fr dataframe with response variable in its first column
#' @param reTrms list with random effect terms
#' @return objective function which maps the covariance parameters \code{theta} to the corresponding (restricted) deviance
#' @export
mkLmerDevfun_R <- function(fr, X, reTrms, REML = TRUE, start = NULL, verbose = 0, ...){
  stopifnot( is.data.frame(fr), NROW(fr) > 0L )
  stopifnot( is.matrix(X), NROW(X) > 0L )

  y <- stats::model.response(fr)   ##fr[, 1L]
  stopifnot( is.numeric(y) )


  Zt <- reTrms[["Zt"]]
  Lambdat <- reTrms[["Lambdat"]]

  theta <- reTrms[["theta"]]
  mapping <- function(theta) theta[reTrms[["Lind"]]]

  n <- length(y)
  p <- NCOL(X)
  q <- NROW(Zt)

  stopifnot( NROW(X) == n, NCOL(Zt) == n)

  # raw weights matrix
  w <- stats::model.weights(fr)
  sqrtW <- if (! is.null(w) && is.numeric(w))
    Matrix::Diagonal(n=n, x=sqrt(w)) else Matrix::Diagonal(n=n)

  # # normalize weights (necessary?)
  # W <- W / sum(W)


  # fixed by design
  WX <- sqrtW %*% X
  Wy <- sqrtW %*% y
  ZtW <- Zt %*% sqrtW

  ZtWX <- ZtW %*% WX
  ZtWy <- ZtW %*% Wy
  XtWX <- Matrix::crossprod(WX)
  XtWy <- Matrix::crossprod(WX, Wy)

  # mutable values stored as closure for function below (realized as local environment)
  local({
    b <- cu <- numeric(q)
    beta <- numeric(p)
    # initial Cholesky factor L_theta
    L <- Matrix::Cholesky(Matrix::tcrossprod(Lambdat %*% ZtW), LDL = FALSE, Imult = 1L)

    Lambdat <- Lambdat                # stored here b/c x slot will be updated
    mu <- numeric(n)                  # conditional mean of response
    RZX <- matrix(0L, nrow=q, ncol=p)  # intermediate matrix in solution
    u <- numeric(q)                   # conditional mode of spherical random effects

    # in order to use lme4::optimizeLmer() we needed some
    # objects like pp in the environment.
    #
    # pp <- do.call(merPredD$new, c(reTrms[c("Zt", "theta",
    #                                        "Lambdat", "Lind")], n = nrow(X), list(X = X)))
    #
    # I follow pureR here; I leave lme4-path at the optimization step

    # profiled (restricted) deviance function
    function(theta) {
      # Step 1: update rel. covariance factor ----
      Lambdat@x[] <<- mapping(theta)
      L <<- Matrix::update(L, Lambdat %*% ZtW, mult = 1L)

      # Step 2: solve normal equation ----
      cu[] <<- as.vector(Matrix::solve(L, Matrix::solve(L, Lambdat %*% ZtWy, system = "P"),
                              system = "L"))
      RZX[] <<- as.vector(Matrix::solve(L, Matrix::solve(L, Lambdat %*% ZtWX, system = "P"),
                                system = "L"))
      RXtRX <- as(XtWX - crossprod(RZX), "dpoMatrix")
      # solutions
      beta[] <<- as.vector(Matrix::solve(RXtRX, XtWy - crossprod(RZX, cu)))
      u[] <<- as.vector(Matrix::solve(L, Matrix::solve(L, cu - RZX %*% beta, system = "Lt"),
                              system = "Pt"))


      # Step 3: update linear predictor ----
      # b = mu_{B | Y=yobs}
      b[] <<- as.vector(Matrix::crossprod(Lambdat, u))
      mu[] <<- as.vector(X %*% beta + Matrix::crossprod(Zt, b))  # + offset ##ZZZ use offset here! where does it come frome?
      wtres <- sqrtW * (y - mu)


      # Step 4: profiled deviance -----
      degFree <- n
      pwrss <- sum(wtres^2) + sum(u^2)
      logDet <- 2L * Matrix::determinant(L, logarithm = TRUE)$modulus

      if (isTRUE(REML)){
        logDet <- logDet + Matrix::determinant(RXtRX, logarithm = TRUE)$modulus
        degFree <- degFree - p
      }
      attributes(logDet) <- NULL

      ## debugging
      if (verbose > 0L){
        cat("theta: ", round(theta, 2), " ")
        cat("beta: ", paste(round(beta, 2), collapse = ", "), "\t")
        cat("b: ", paste(round(b, 2), collapse = ", "), "\t")
        cat("s2: ", round(pwrss/degFree,2), "\n")
      }

      # profiled deviance
      objVal <- logDet + degFree * (1L + log(2 * pi * pwrss) - log(degFree))

      ## add corresponding model parameter fits as attributes
      attr(objVal, which = "beta") <- beta
      attr(objVal, which = "b") <- b
      attr(objVal, which = "resVar") <- pwrss/degFree

      objVal
    }
  })

}



