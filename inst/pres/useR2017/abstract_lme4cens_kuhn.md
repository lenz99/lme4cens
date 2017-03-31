---
title: "Censoring in random effects models"
author: |
  | Matthias Kuhn^1^ and Ingo Roeder^1,2^
  |
  | 1. Institute for Medical Informatics and Biometry (IMB), Faculty of Medicine Carl Gustav Carus
  | Technische Universität Dresden, Dresden, Germany
  |
  | 2. National Center for Tumor Diseases, Partner Site Dresden, Germany
output:
  html_document: default
  pdf_document: default
  word_document: default
references:
- URL: https://www.jstatsoft.org/index.php/jss/article/view/v067i01
  author:
  - family: Bates
    given: Douglas
  - family: Mächler
    given: Martin
  - family: Bolker
    given: Ben
  - family: Walker
    given: Steve
  container-title: Journal of Statistical Software
  id: lmmLme4Bates2015
  issue: 1
  issued:
    month: 10
    year: 2015
  page: 1-48
  title: Fitting Linear Mixed-Effects Models Using lme4
  volume: 67
- URL: http://cran.r-project.org/package=censReg
  author:
  - family: Henningsen
    given: Arne
  id: henningsencensreg
  issued:
    year: 2010
  title: 'censReg: Censored Regression (Tobit) Models. R Package'
  type: article

---

**Keywords**: censoring, random effects model

**Webpages**: https://github.com/lenz99-/lme4cens


Random effects models are one option in regression analysis when the data has a hierarchical structure, i.e. when the observations are _not_ independent and identically distributed.
Random effects models account for the induced correlation by assigning zero-centered normally distributed random variables to different hierarchy levels in the data and by predicting their value during the fitting process.
As an example, repeated measurements in a longitudinal study are nested within subjects and e.g. a simple random effects model would provide for a random intercept per subject.
In the *R* world, random effects models are handled -- among others -- in the **lme4** package.
Given a data sample, the parameters of a random effects model are estimated in **lme4** via (restricted) maximum likelihood.

If a response is not known exactly but only to have occurred within a certain interval we speak of censoring.
It is a well-known concept from time-to-event analysis but censoring also occurs e.g. when the response is a concentration and the measuring device has a lower detection limit.
There exist *R*-implementations for regression analysis with censored response, most notably the **survival**-package with its `survreg`-function for parametric regression models.
Regarding random effects models, there are custom *R*-packages (e.g. **censReg**, @henningsencensreg) that support censoring as well, at least for simple random intercept models.

The **lme4**-package is modular and developers are encouraged to reuse functionality for model enhancements or specializations [@lmmLme4Bates2015].
We have developed our *R*-package **lme4cens** that builds on **lme4** to fit simple random effects models with a censored response.
In particular, the re-use of the formula-module facilitates model specification.
The censoring information is encoded via **survival**'s `Surv`-object that allows for a flexible specification of (a combination of) left-, right- and interval-censored responses with varying censoring levels.


# References
