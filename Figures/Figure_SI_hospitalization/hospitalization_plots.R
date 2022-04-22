
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
library(wesanderson)

options(width=100)
options(scipen=100)

enableJIT(3) 

source("functions.R")

set.seed(101)

## Read in parameters and set initial conditions
exog_parms <- read.csv("exog_parms_main.csv")
exog_parms$final.time <- 548
final.time = exog_parms$final.time

discount_rate = 0.04                                              # annual discount rate
discount_factor = (1/(1+discount_rate))^(1/365)       # daily discount factor

phi <- 0.8554632
model_prod_loss <- ((1-phi)*100)

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

## Read in lockdowns
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

## Read in decentralized value and policy functions
solved_values = read.csv("../../Results/Figure2_Models_epi/Data_For_plots/value_functions/nextgen_eqm_vpfn.csv")

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

## Generate time path of hospitalization under baseline policy

setwd("../../Results/Figure2_Models_epi/Data_For_plots/time_series")
main_eqm_path <- read_csv("planner_35pt.csv")

# Baseline hospitalization
hospitalized <- rep(NA,length.out=nrow(main_eqm_path))
hospitalized[1] <- 0
for(i in 1:(nrow(main_eqm_path)-1)) {
	hospitalized[i+1] <- hospitalized[i]*(P_HH - P_HR - P_HD) + main_eqm_path$I[i]*P_IH 
}

main_eqm_path <- cbind(main_eqm_path, hospitalized, type="targeted_isolation")

# Lockdown hospitalization
hospitalized <- rep(NA,length.out=nrow(results_policy))
hospitalized[1] <- 0
for(i in 1:(nrow(results_policy)-1)) {
	hospitalized[i+1] <- hospitalized[i]*(P_HH - P_HR - P_HD) + results_policy$I[i]*P_IH 
}

results_policy <- cbind(results_policy, hospitalized, type="lockdown")

hosp_series_small <- rbind(main_eqm_path, results_policy)
hosp_series_small <- hosp_series_small[,c("time","hospitalized","type")]


implied_hosp_rate <- ggplot(data = hosp_series_small, aes(x = time, group = type, linetype = type)) + 
	geom_line(aes(y = hospitalized), size = 1.1) + 
	theme_bw() +
	labs(y = "% of population", x = "Day", title = "Implied hospitalization rate") +
	scale_linetype_manual("", values = c('targeted_isolation' = 'solid', 'lockdown' = 'dotted'), labels = c("\nTargeted\nisolation", "\nBlanket\nlockdown\n"))

# hospital_plot

png(paste0("../../../Figure_SI_hospitalization/hospitalization2.png"), width=700, height=600)
# hospital_plot
implied_hosp_rate
dev.off()





# main_eqm_path_small <- pivot_longer(main_eqm_path[,c("time","I","hospitalized")], cols = c(I, hospitalized), names_to="outcome", values_to="value")

# implied_hosp_rate <- ggplot(data = main_eqm_path_small, aes(x = time, group = outcome, linetype = outcome)) + 
# 	geom_line(aes(y = value), size = 0.75) + 
# 	theme_bw() +
# 	labs(y = "% of population", x = "Day", title = "Implied hospitalization rate") +
# 	scale_linetype_manual("", values = c('hospitalized' = 'solid', 'I' = 'dashed'), labels = c("\nHospitalization\nrate", "\nInfection\nprevalence\n"))













##########################################################################################
# OLD HEATMAP STUFF
##########################################################################################



## Read in rhoc data
# setwd("../../../Figure3_Sensitivity_analysis/data/ngm_11_rhoc")

# ts_filenames = list.files(pattern="*.csv")

# time_series_list = list()
# for(name in seq_along(ts_filenames)) {
# 	label = gsub(".csv", "", ts_filenames[name])
# 	label = gsub("_0_0.+", "", label)
# 	rhoc = str_extract(label, "(?<=rhoc_)\\d+")
# 	phi = str_extract(label, "(?<=phi_)\\d+[:punct:]\\d+")

# 	label = gsub("_rhoc_+.+", "", label)
# 	time_series_list[[name]] = cbind(read_csv(paste0(ts_filenames[name])),
# 	                      type=paste0(label), rhoc=as.numeric(rhoc), phi=as.numeric(phi))
# 	## Generate hospitalization series for all paths
# 	hospitalized <- rep(NA,length.out=nrow(time_series_list[[name]]))
# 	hospitalized[1] <- 0
# 	for(i in 1:(nrow(time_series_list[[name]])-1)) {
# 		hospitalized[i+1] <- hospitalized[i]*(P_HH - P_HR - P_HD) + time_series_list[[name]]$I[i]*P_IH 
# 	}

# 	time_series_list[[name]] = cbind(time_series_list[[name]], hospitalized=hospitalized)
# }

# ts_dfrm = rbindlist(time_series_list)

# scenarios_summary = ts_dfrm %>% 
# 	group_by(type,rhoc,phi) %>%
# 	filter(type == "ngm_eqm") %>%
# 	summarise(max_hospital_burden = max(hospitalized)*100)

# cl_data =  pivot_wider(scenarios_summary, id_cols=c("type","rhoc","phi"), names_from=type, values_from=c(max_hospital_burden), names_sep="_") %>%
# 	mutate(productivity_loss = (1 - phi)*100) %>%
# 	mutate(c_l_contact_ratio = rhoc/7.513051)

# cl_model_ratio <- 5.166426/7.513051

# ## Read in rhoo data
# setwd("../ngm_10_rhoo")

# ts_filenames = list.files(pattern="*.csv")

# time_series_list = list()
# for(name in seq_along(ts_filenames)) {
# 	label = gsub(".csv", "", ts_filenames[name])
# 	label = gsub("_0_0.+", "", label)
# 	rhoo = str_extract(label, "(?<=rhoo_)\\d+")
# 	phi = str_extract(label, "(?<=phi_)\\d+[:punct:]\\d+")

# 	label = gsub("_rhoo_+.+", "", label)
# 	time_series_list[[name]] = cbind(read_csv(paste0(ts_filenames[name])),
# 	                      type=paste0(label), rhoo=as.numeric(rhoo), phi=as.numeric(phi))
# 	## Generate hospitalization series for all paths
# 	hospitalized <- rep(NA,length.out=nrow(time_series_list[[name]]))
# 	hospitalized[1] <- 0
# 	for(i in 1:(nrow(time_series_list[[name]])-1)) {
# 		hospitalized[i+1] <- hospitalized[i]*(P_HH - P_HR - P_HD) + time_series_list[[name]]$I[i]*P_IH 
# 	}

# 	time_series_list[[name]] = cbind(time_series_list[[name]], hospitalized=hospitalized)
# }

# ts_dfrm = rbindlist(time_series_list)

# scenarios_summary = ts_dfrm %>% 
# 	group_by(type,rhoo,phi) %>%
# 	filter(type == "ngm_eqm") %>%
# 	summarise(max_hospital_burden = max(hospitalized)*100)

# ua_data =  pivot_wider(scenarios_summary, id_cols=c("type","rhoo","phi"), names_from=type, values_from=c(max_hospital_burden), names_sep="_") %>%
# 	mutate(productivity_loss = (1 - phi)*100) %>%
#     mutate(unavoid_avoid_contact_ratio = rhoo/(5.166426+7.513051) )

# ua_model_ratio <- 3.548544/(5.166426+7.513051)

# ## Compute measures for heatmap

# # cl_hospital = ggplot(data = cl_data, aes(x=productivity_loss, y=c_l_contact_ratio)) +
# #             geom_tile(aes(fill = ngm_eqm)) +
# #             theme_classic() +
# #             labs(y = "Ratio of consumption to labor contacts", x = "Productivity loss from infection (%)", title = "Maximum hospital burden", fill = "% of pop")

#  cl_hospital <- ggplot(cl_data) + 
# 		geom_contour_filled(aes(x=productivity_loss,y=c_l_contact_ratio, z = ngm_eqm), size=0.009) + 
# 		scale_fill_manual(values = wes_palette(14, name = "Zissou1", type = "continuous"), name = "% of pop") + 
# 		theme_minimal() + 
# 		ylab("Ratio of consumption to labor contacts") + xlab("Productivity loss from infection (%)") + ggtitle("Maximum hospital burden") + 
# 		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_text(color="grey30")) + 
# 		#scale_y_continuous(breaks = seq(0,9,by = 1)) +
# 		#scale_x_continuous(breaks = seq(0,90,by = 10)) +
# 		geom_point(aes(x=model_prod_loss, y=cl_model_ratio), fill="white", colour="#4f418b", shape=21,size = 1.5, stroke =1)

# # ua_hospital = ggplot(data = ua_data, aes(x=productivity_loss, y=unavoid_avoid_contact_ratio)) +
# #             geom_tile(aes(fill = ngm_eqm)) +
# #             theme_classic() +
# #             labs(y = "Ratio of unavoidable to avoidable contacts", x = "Productivity loss from infection (%)", title = "Maximum hospital burden", fill = "% of pop")

# ua_hospital <- ggplot(ua_data) + 
# 		geom_contour_filled(aes(x=productivity_loss, y=unavoid_avoid_contact_ratio, z = ngm_eqm), size=0.009) + 
# 		scale_fill_manual(values = wes_palette(14, name = "Zissou1", type = "continuous"), name = "% of pop") + 
# 		theme_minimal() + 
# 		ylab("Ratio of unavoidable to avoidable contacts") + xlab("Productivity loss from infection (%)") + ggtitle("Maximum hospital burden") + 
# 		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_text(color="grey30")) + 
# 		#scale_y_continuous(breaks = seq(0,9,by = 1)) +
# 		#scale_x_continuous(breaks = seq(0,90,by = 10)) +
# 		geom_point(aes(x=model_prod_loss, y=ua_model_ratio), fill="white", colour="#4f418b", shape=21,size = 1.5, stroke =1)

# ## Make plots

# hospital_plot <- implied_hosp_rate / (ua_hospital | cl_hospital) + plot_annotation(tag_levels = "a")
