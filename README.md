# lme4cens: Simple Linear Mixed Effect Models and Censoring

The *R*-package **lme4cens** builds on **lme4** to fit simple random effects models with a censored response.
It re-uses the formula-module from **lme4** to facilitate model specification.
The censoring information is encoded via **survival**'s `Surv`-object that allows for a flexible specification of
(a combination of) left-, right- and interval-censored responses with flexible censoring levels per observation.

The random effect structure is currently limited to the most simple case, namely models with a single random intercept.
Model fitting is via maximum likelihood (`ML`), residual maximum likelihood (`REML`) is not supported.

A good choice of starting values is helpful, although there is a heuristic in place if none are given.
The fitted parameter values may depend on the choice of starting values.
As with all non-trivial optimization problems it is good practice to check convergence with different start values.
