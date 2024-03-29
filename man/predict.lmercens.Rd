% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{predict.lmercens}
\alias{predict.lmercens}
\title{Predict method for \code{lmercens} objects.}
\usage{
\method{predict}{lmercens}(object, newdata = NULL, re.form = NULL, ...)
}
\arguments{
\item{object}{\code{lmercens} model object}

\item{newdata}{dataframe with covariate values for prediction. Default is \code{NULL} which is to fall back to the training data.}

\item{re.form}{right-side formula for random effects to condition on in case of prediction on training data. If \code{NULL}, include all random effects; if \code{~0} or \code{NA}, include no random effects.}
}
\description{
Predict method for \code{lmercens} objects.
}
