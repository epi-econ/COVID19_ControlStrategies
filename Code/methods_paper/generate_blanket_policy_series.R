# Script to make some time series outputs of imposing a blanket policy. If this works, let it take the ploicy parameters as an input, so that the script can be embedded into figure2_lines_sketch.R

library(ggplot2)
library(patchwork)
library(doSNOW)
library(doParallel)
library(optimParallel)
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

parms <- read.csv("../../Results/Figure2_Models_epi/Data_For_plots/value_functions/calibrated_parameters.csv")
parms$kappa_c <- 1
parms$kappa_l <- 1

discount_rate = 0.04                                              # annual discount rate
discount_factor = (1/(1+discount_rate))^(1/365)       # daily discount factor

## Read in summary of losses
policy_summary = read.csv("../../Results/Figure_SI_blanket_sensitivity/ensemble_summary.csv")
best_policy_in_disease_terms = which.min(policy_summary$cases_per_100k[(1:(nrow(policy_summary)-2))])
best_policy_in_economic_terms = which.min(policy_summary$PV_losses[(1:(nrow(policy_summary)-2))])

## Set blanket policy parameters
policy_list = list()
policy_list$switch = 1 #1 if policy will be put in place
policy_list$start_date = policy_summary$start_date[best_policy_in_disease_terms] #policy start time in days
policy_list$end_date = policy_summary$end_date[best_policy_in_disease_terms] #policy end time in days (2-week shutdown)
policy_list$labor_supply_level = policy_summary$labor_supply_level[best_policy_in_disease_terms] #reduce labor supply to this much of initial (day 1) level
policy_list$start_trigger = 0.1 #start 1-day stay-at-home if infection crosses this threshold

best_params = as.data.frame(unlist(policy_list),policy_summary$cases_per_100k[best_policy_in_disease_terms],policy_summary$PV_losses[best_policy_in_disease_terms])

write.csv(best_params, file = "~/Dropbox/Corona_epicon/Results/Figures_data_FINAL/best_disease_policy_parameters.csv")

## Read in decentralized value and policy functions
solved_values = read.csv("../../Results/Figure2_Models_epi/Data_For_plots/value_functions/nextgen_eqm_vpfn.csv")

## Read in parameters and set initial conditions
exog_parms = read.csv("exog_parms_main.csv")
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

## Generate blanket policy time series
results_policy = generate_time_series_policy(solved_values,policy_list)

## Save output
write_csv(results_policy, path = "../../Results/Figure2_Models_epi/Data_For_plots/blanket_1.csv")

## Generate no policy, append policy indicator to results_policy
nopolicy <- cbind(generate_time_series(solved_values=solved_values),policy="no")
results_policy = cbind(results_policy,policy="yes")
results_planner <- cbind(read_csv("../../Results/Figure2_Models_epi/Data_For_plots/planner_35pt.csv"),policy="planner")

## Append some columns with contacts stuff

results_policy <- results_policy %>%
	mutate(SI_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_I, l_S = labor_S, l_I = labor_I) )) %>%
	mutate(SR_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_R, l_S = labor_S, l_I = labor_R) )) %>%
	mutate(RI_contacts = contacts(parms, choices = data.frame(c_S = consumption_R, c_I = consumption_I, l_S = labor_R, l_I = labor_I) )) %>%
	mutate(SS_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_S, l_S = labor_S, l_I = labor_S) )) %>%
	mutate(RR_contacts = contacts(parms, choices = data.frame(c_S = consumption_R, c_I = consumption_R, l_S = labor_R, l_I = labor_R) )) %>%
	mutate(prob_contact_I = (SI_contacts + RI_contacts)/(SI_contacts + RI_contacts + SR_contacts + SS_contacts + RR_contacts)) %>%
	mutate(prob_contact_I_weighted = (SI_contacts*S*I + RI_contacts*R*I)/(SI_contacts*S*I + RI_contacts*R*I + SR_contacts*S*R + SS_contacts*S*S + RR_contacts*R*R))

nopolicy <- nopolicy %>%
	mutate(SI_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_I, l_S = labor_S, l_I = labor_I) )) %>%
	mutate(SR_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_R, l_S = labor_S, l_I = labor_R) )) %>%
	mutate(RI_contacts = contacts(parms, choices = data.frame(c_S = consumption_R, c_I = consumption_I, l_S = labor_R, l_I = labor_I) )) %>%
	mutate(SS_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_S, l_S = labor_S, l_I = labor_S) )) %>%
	mutate(RR_contacts = contacts(parms, choices = data.frame(c_S = consumption_R, c_I = consumption_R, l_S = labor_R, l_I = labor_R) )) %>%
	mutate(prob_contact_I = (SI_contacts + RI_contacts)/(SI_contacts + RI_contacts + SR_contacts + SS_contacts + RR_contacts)) %>%
	mutate(prob_contact_I_weighted = (SI_contacts*S*I + RI_contacts*R*I)/(SI_contacts*S*I + RI_contacts*R*I + SR_contacts*S*R + SS_contacts*S*S + RR_contacts*R*R))

results_planner <- results_planner %>%
	mutate(SI_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_I, l_S = labor_S, l_I = labor_I) )) %>%
	mutate(SR_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_R, l_S = labor_S, l_I = labor_R) )) %>%
	mutate(RI_contacts = contacts(parms, choices = data.frame(c_S = consumption_R, c_I = consumption_I, l_S = labor_R, l_I = labor_I) )) %>%
	mutate(SS_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_S, l_S = labor_S, l_I = labor_S) )) %>%
	mutate(RR_contacts = contacts(parms, choices = data.frame(c_S = consumption_R, c_I = consumption_R, l_S = labor_R, l_I = labor_R) )) %>%
	mutate(prob_contact_I = (SI_contacts + RI_contacts)/(SI_contacts + RI_contacts + SR_contacts + SS_contacts + RR_contacts)) %>%
	mutate(prob_contact_I_weighted = (SI_contacts*S*I + RI_contacts*R*I)/(SI_contacts*S*I + RI_contacts*R*I + SR_contacts*S*R + SS_contacts*S*S + RR_contacts*R*R))

results <- rbind(results_policy,nopolicy,results_planner)

## Diagnostic plots

plot_base <- ggplot(data=results, aes(x=time, group=policy, color=policy))

infection = plot_base + 
            geom_line(aes(y=I*100,linetype=policy), size=1) +
            theme_bw() + ggtitle("Infecteds") +
            xlab("Day") + ylab("Proportion of population (%)") +
            geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
            geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
            scale_size_discrete(range = c(1,1.5)) +
            guides(linetype=FALSE, size=FALSE)

infection_zoomed = ggplot(data=(results %>% filter(time<(as.numeric(best_params["end_date",])+28), time>(as.numeric(best_params["start_date",])-14))), aes(x=time)) +
            geom_line(aes(y=I*100,linetype=policy,color=policy), size=1) +
            theme_bw() + ggtitle("Infecteds: Day 45-75") +
            xlab("Day") + ylab("Proportion of population (%)") + 
            geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
            geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
            scale_size_discrete(range = c(1,1.5)) +
            guides(linetype=FALSE, size=FALSE)

total_contacts = plot_base + 
            geom_line(aes(y=consumption_contacts+labor_contacts+other_contacts), size=1) +
            geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
            geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
            theme_bw() + ggtitle("Total S-I contacts") +
            xlab("Day") + ylab("Daily contacts") + 
            guides(color=FALSE)

Reff = plot_base + 
            geom_line(aes(y=Reff,linetype=policy), size=1) +
            theme_bw() + ggtitle("R effective") +
            xlab("Day") + ylab("R_t") + 
            geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
            geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
            scale_size_discrete(range = c(1,1.5)) +
            guides(linetype=FALSE, size=FALSE)

prob_contact_wt = plot_base + 
            geom_line(aes(y=prob_contact_I_weighted,linetype=policy), size=1) +
            theme_bw() + ggtitle("Probability random S or R contacts an I") +
            xlab("Day") + ylab("Probability") + 
            geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
            geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
            scale_size_discrete(range = c(1,1.5)) +
            guides(linetype=FALSE, size=FALSE)

weighted_labor_supplies = plot_base +
            geom_line(aes(y = S*labor_S, group=policy, linetype=policy), size=1, color = "darkgreen") +
            geom_line(aes(y = I*labor_I, group=policy, linetype=policy), size=1, color = "firebrick4") +
            geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
            geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
            theme_bw() + ggtitle("Aggregate labor supplies\n(green = S, red = I)") +
            xlab("Day") + ylab("Normalized person-hours")

blanket_plot = ((infection | total_contacts | Reff) / (prob_contact_wt |infection_zoomed | weighted_labor_supplies)) + 
      plot_layout(heights = c(2,2)) + 
      plot_annotation(
        title = 'Diagnostic: assessing the example general lockdown',
        subtitle = "policy = lockdown"
      )

blanket_plot

png(paste0("../../Results/Figures_data_FINAL/general_lockdown_plot.png"), width=700, height=600)
blanket_plot
dev.off()
