
#' Sleepstudy data with additional censoring information (event2).
#' @source lme4 package
"sleepstudy"

#' Balanced subset of Sleepstudy with only four patients, including censoring information (event2).
#'
#' @details
#' For testing purposes with `interval`-type censoring status
#'
#' * `event2` [numeric] status variable
#' * `event3` [numeric] status variable with clear-cut censoring boundaries (<312 L and > 350 R)
#'
#' @source lme4 package
"sleepstudy2"

#' Fair's Extramartial Affairs Data.
#'
#' Cross-section data from a survey conducted by Psychology Today in 1969.
#' An event-column was added for easier use as [survival::Surv]-object.
#'
#' @details
#' The variables in the dataframe are
#'
#' * `affairs` [numeric] How often engaged in affairs during the past year?
#' * `gender` [factor] indicating gender.
#' * `age` [numeric] coding age in years.
#' * `religiousness` [numeric] coding religiousness from 1=anti, to 5=very
#' * ...
#' * `event` [numeric] status-coding, 0: left censored, 1: event
#' * `event2` [numeric] status coding for left- and right-censoring, 0: right-censored, 1: event, 2: left-censored
#'
#' @source AER package
"Affairs"
