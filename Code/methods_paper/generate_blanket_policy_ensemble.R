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
library(stargazer)
library(tidyverse)

options(width=100)
options(scipen=100)

enableJIT(3) 

source("functions.R")

## Read in SIR path
SIR = read.csv("../../Results/Figure_SI_blanket_sensitivity/nextgen_SIR.csv")

## Read in lockdowns
lockdowns = read.csv("../../Results/Figure_SI_blanket_sensitivity/lockdowns.csv")

set.seed(101)

ncores = 16
ensemble_length = 2500

## Distribution of start dates
max_date = 100
start_date_vec = seq(1:max_date)
start_date_probs = SIR$I[1:max_date]
start_date_probs = start_date_probs/sum(start_date_probs)

## Distribution of lockdown lengths
# lockdown_lengths = pmin(ceiling(rexp(n=ensemble_length, rate=1/21)),60)
#lockdown_lengths = base::sample(lockdowns$duration_till_first_reopen, size=ensemble_length, replace=TRUE)
lockdown_lengths = round(rtrunc(ensemble_length, "exp", a=min(lockdowns$duration_till_first_reopen), b=max(lockdowns$duration_till_first_reopen), rate=1/mean(lockdowns$duration_till_first_reopen)),0)
lockdown_lengths
summary(lockdown_lengths)

## Draw blanket policy parameters
start_dates = sample(start_date_vec, size=ensemble_length, prob=start_date_probs, replace=TRUE) #c(30,45,60) # policy start time in days since epidemic start
end_dates = start_dates + lockdown_lengths # policy end time in days (14 = 2-week shutdown)
labor_supply_level = runif(ensemble_length, 0, 1) #rbeta(n=ensemble_length, shape1=3, shape2=1.5) # reduce labor supply to this much of initial (day 1) level

blanket_lockdowns = data.frame(start=start_dates, end=end_dates, severity=labor_supply_level)
blanket_lockdowns
summary(blanket_lockdowns)

summary_stats_for_paper <- list(empirical_duration = lockdowns$duration_till_first_reopen, empirical_timing = lockdowns$duration_till_first_close, simulated_duration = (blanket_lockdowns$end - blanket_lockdowns$start), simulated_timing = blanket_lockdowns$start, simulated_severity = blanket_lockdowns$severity)
attributes(summary_stats_for_paper) <- list(names = names(summary_stats_for_paper),
    row.names=1:2500, class='data.frame')
summary(summary_stats_for_paper)

stargazer(summary_stats_for_paper)

ggplot(data = blanket_lockdowns) + geom_jitter(aes(x=start_dates,y=end_dates,color=severity), width=0.1, height=0.1)

## Create blanket policy list object
policy_list_ensemble = list()
policy_dfrm_ensemble = as.data.frame(matrix(nrow=ensemble_length,ncol=4))
colnames(policy_dfrm_ensemble) = c("start_date","end_date","labor_supply_level","ensemble_id")
for(i in 1:ensemble_length) {
	policy_list = list()
	policy_list$switch = 1 #1 if policy will be put in place
	policy_list$start_date = start_dates[i] 
	policy_list$end_date = end_dates[i]  
	policy_list$labor_supply_level = labor_supply_level[i]
	policy_list_ensemble[[i]] = policy_list
	policy_dfrm_ensemble[i,] = c(unlist(policy_list)[-1],i)
}

## Read in decentralized value and policy functions
solved_values = read.csv("../../Results/Figure2_Models_epi/Data_For_plots/value_functions/nextgen_eqm_vpfn.csv")

## Read in parameters and set initial conditions
exog_parms <- read.csv("../../Results/Figure2_Models_epi/Data_For_plots/value_functions/calibrated_parameters.csv")
exog_parms$final.time <- 548
final.time = exog_parms$final.time

discount_rate = 0.04                                              # annual discount rate
discount_factor = (1/(1+discount_rate))^(1/365)       # daily discount factor

total_population = 331002651				# 331002651 is the US population
I_0 = 331.002651/total_population			# set to get 1-in-a-million as an initial condition
S_0 = 1 - I_0
R_0 = 0
SIR_init = data.frame(S=S_0, I=I_0, R=R_0)

## Create grid and other necessary objects
S_gridlength = length(unique(solved_values$S))
I_gridlength = length(unique(solved_values$I))
R_gridlength = length(unique(solved_values$R))
grid_pieces = build_grid(S_gridlength, I_gridlength, R_gridlength)
grid_list = grid_pieces[1:3]

grid_dfrm = grid_pieces$dfrm
simplex_check = rep(1,length.out=nrow(grid_dfrm))

## Generate blanket policy time series ensemble
results_policy_ensemble = list()
# results_nopolicy_ensemble = list()
# for(i in 1:ensemble_length){
# 	results_policy_ensemble[[i]] = cbind(generate_time_series_policy(solved_values=solved_values,policy_list=policy_list_ensemble[[i]]),ensemble_id=i)
# 	results_nopolicy_ensemble[[i]] = cbind(generate_time_series(solved_values=solved_values),ensemble_id=i)
# 	message("Finished virtual copy ",i)
# }

# policy_ensemble_dfrm = rbindlist(results_policy_ensemble)
# nopolicy_ensemble_dfrm = rbindlist(results_nopolicy_ensemble)

cl = makeCluster(ncores)
registerDoParallel(cl)

policy_ensemble_dfrm <- foreach(i=1:ensemble_length, .export=ls(globalenv()), .packages=c("chebpol"), .combine=rbind, .inorder=TRUE) %dopar% {
	virtual_copy = cbind(generate_time_series_policy(solved_values=solved_values,policy_list=policy_list_ensemble[[i]]),ensemble_id=i)
	virtual_copy
}

stopCluster(cl)

# View(policy_ensemble_dfrm)

## Diagnostic plots
# bigthing <- rbind(cbind(policy_ensemble_dfrm,policy="yes"), cbind(nopolicy_ensemble_dfrm,policy="no"))

# plot_base <- ggplot(data=bigthing, aes(x=time,linetype=policy,color=policy))

# infection <- plot_base + 
# 	geom_line(aes(y=I*100), size=1) +
# 	geom_vline(xintercept=policy_list$start_date, linetype="dashed", alpha=0.5) +
# 	geom_vline(xintercept=policy_list$end_date, linetype="dashed", alpha=0.5) +
#     theme_bw() + ggtitle("Infecteds") +
#     xlab("Day") + ylab("Proportion of population (%)") + 
#     guides(linetype=FALSE, size=FALSE)

# total_contacts = plot_base + 
#     geom_line(aes(y=consumption_contacts+labor_contacts+other_contacts), size=1) +
# 	geom_vline(xintercept=policy_list$start_date, linetype="dashed", alpha=0.5) +
# 	geom_vline(xintercept=policy_list$end_date, linetype="dashed", alpha=0.5) +
#     theme_bw() + ggtitle("Total S-I contacts") +
#     xlab("Day") + ylab("Daily contacts") + 
#     guides(color=FALSE)

# infection / total_contacts

## Save output
write_csv(policy_ensemble_dfrm, path = "../../Results/Figure_SI_blanket_sensitivity/data/ensemble.csv")
write_csv(policy_dfrm_ensemble, path = "../../Results/Figure_SI_blanket_sensitivity/ensemble_policy_parameters.csv")
