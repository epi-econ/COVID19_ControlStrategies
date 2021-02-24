# Script to make some time series outputs of imposing a blanket policy. If this works, let it take the ploicy parameters as an input, so that the script can be embedded into figure2_lines_sketch.R

library(ggplot2)
library(patchwork)
library(doSNOW)
library(doParallel)
library(optimParallel)
library(truncdist)
library(chebpol)
library(pracma)
library(rootSolve)
library(compiler)
library(data.table)
library(dplyr)
library(fields)
library(progress)
library(tidyverse)

options(width=100)
options(scipen=100)

enableJIT(3) 

source("functions.R")

set.seed(101)
## Read in decentralized and planner time paths, and decentralized policy function
eqm_path = read.csv("../../Results/Figure2_Models_epi/Data_For_plots/eqm_35pt.csv")
plan_path = read.csv("../../Results/Figure2_Models_epi/Data_For_plots/planner_35pt.csv")
# eqm_vpfn = read.csv("../../Results/Figure2_Models_epi/Data_For_plots/value_functions/eqm_vpfn.csv")

## Read in parameters and set initial conditions
exog_parms <- read.csv("exog_parms_main.csv")
exog_parms$final.time <- 548
final.time = exog_parms$final.time

discount_rate = 0.04                                              # annual discount rate
discount_factor = (1/(1+discount_rate))^(1/365)       # daily discount factor

total_population = 331002651				# 331002651 is the US population
I_0 = 331.002651/total_population			# set to get 1-in-a-million as an initial condition
S_0 = 1 - I_0
R_0 = 0
SIR_init = data.frame(S=S_0, I=I_0, R=R_0)

## Set hospitalization parameters
P_IH <- 0.199		 					# 19.9% hospitalized, per our calculations in SI
P_R <- exog_parms$pi_r
P_D <- exog_parms$pi_d
P_I <- 1 - P_R - P_D
hospital_stay_duration <- 10 			# 10 days, from Wang et al (2020)
P_HH <- 1/hospital_stay_duration 		# equal probability of exiting hospital each day
est_cfr <- 0.015
P_HR <- (1-P_HH)*(1-est_cfr)
P_HD <- (1-P_HH)*est_cfr

## Calculate parameters to make things consistent
P_IR <- P_R - P_IH*P_HR
P_ID <- P_D - P_IH*P_HD
P_II <- P_I - P_IH*P_HH

P_II
P_IR
P_ID

## Check adding up
P_II + P_IH + P_IR + P_ID
P_HH + P_HR + P_HD

## Generate hospitalization series for paths under decentralized behavior
hospitalized <- rep(NA,length.out=nrow(eqm_path))
hospitalized[1] <- 0
for(i in 1:(nrow(eqm_path)-1)) {
	hospitalized[i+1] <- hospitalized[i]*(P_HH - P_HR - P_HD) + eqm_path$I[i]*P_IH 
}

eqm_path <- data.frame(eqm_path,hospitalized = hospitalized)

## Plot path
hospitalization <- ggplot(data = eqm_path, aes(x=time)) + 
	geom_line(aes(y=hospitalized*100), size = 1) +
	labs(y = "% of population", x = "Day", title = "Implied hospitalization rate") + 
	theme_bw()
hospitalization_context <- ggplot(data = eqm_path, aes(x=time)) + 
	geom_line(aes(y=hospitalized*100), size = 1) +
	geom_line(aes(y=I*100), size = 1, linetype="dashed") +
	labs(y = "% of population", x = "Day", title = "Implied hospitalization rate (solid)\ngiven infection rate (dashed)") + 
	theme_bw()

hospitalization_plot <- (hospitalization | hospitalization_context)

png(paste0("../../Results/Figure_SI_hospitalization/hospitalization.png"), width=800, height=500)
hospitalization_plot
dev.off()


## Draw alternate hospitalization parameters

