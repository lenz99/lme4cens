
#' Sleepstudy data with additional censoring information (event2).
#' @source lme4 package
"sleepstudy"

#' Subset of `sleepstudy` data with only four patients,
#' including possible censoring information (`event2` and `event3`).
#'
#' @details
#' This subset of data is  provided with different event columns for testing purposes with double censored censoring status.
#'
#' * `event2` [numeric] status variable with __random__ left- and right-censoring (left cens possible below <300, right cens possible >310)
#' * `event3` [numeric] status variable with __deterministic__ clear-cut left- and right-censoring (boundaries are <212 for Left and >350 for Right cens)
#' * `Reaction3` [numeric] response variable with the clear-cut censoring boundaries (from `event3`) as the response for censored observations
#'
#' @source `lme4` package
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
