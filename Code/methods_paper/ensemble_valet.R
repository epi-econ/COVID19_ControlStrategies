## R script to run through the gamut of ensemble stuff: generate the scenarios, sumamrize them and do the summary plot, then pick the best one of the lot and sketch out that time series

setwd("~/Dropbox/Corona_epicon/Code/methods_paper")
source("generate_blanket_policy_ensemble.R", print.eval=TRUE)

rm(list=ls())
setwd("~/Dropbox/Corona_epicon/Code/methods_paper")
source("ensemble_lines_sketch.R", print.eval=TRUE)

rm(list=ls())
setwd("~/Dropbox/Corona_epicon/Code/methods_paper")
source("generate_blanket_policy_series.R", print.eval=TRUE)
