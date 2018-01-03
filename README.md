
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/medrc)](https://cran.r-project.org/package=medrc) [![Build Status](https://travis-ci.org/DoseResponse/medrc.svg?branch=master)](https://travis-ci.org/DoseResponse/medrc) [![Downloads](https://cranlogs.r-pkg.org/badges/medrc)](https://cranlogs.r-pkg.org/)

medrc
=====

Overview
--------

The medrc package is an extension to the package drc that allows you to fit hierarchical dose-response models and provides inferenc for derived parameters like the effective dose or benchmark doses. Parameters are estimated with a two-stage approach using the package metafor, or using the Lindstrom-Bates algorithm implemented in the package nlme.

Installation
------------

``` r
## You can install medrc from GitHub
# install.packages("devtools")
## first installing drc and drcData
devtools::install_github("daniel-gerhard/drc")
devtools::install_github("daniel-gerhard/drcData")
## then installing the development version of medrc
devtools::install_github("daniel-gerhard/medrc")
```
