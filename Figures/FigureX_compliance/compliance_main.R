#### Script with various information frictions

library(ggplot2)
library(patchwork)
library(doSNOW)
library(doParallel)
library(optimParallel)
library(chebpol)
library(pracma)
library(rootSolve)
library(nleqslv)
library(compiler)
library(data.table)
library(dplyr)
library(fields)
library(progress)
library(tidyverse)
library(viridis)
library(janitor)
library(ggpubr)
library(scales)

options(width=100)
options(scipen=100)

enableJIT(3) 

setwd("~/Dropbox/Corona_epicon/Code/methods_paper/")

source("functions.R")

#####
# Set parameters
#####

discount_rate <- 0.04                                              # annual discount rate
discount_factor <- (1/(1+discount_rate))^(1/365)       # daily discount factor

total_population <- 331002651
I_0 = 331.002651/total_population               # set to get 1-in-a-million as an initial condition
S_0 = 1 - I_0
R_0 = 0
SIR_init = data.frame(S=S_0, I=I_0, R=R_0)
maximum_gain <- 0.91
max_averted_cases <- 1 - 61629/61826

final.time <- 548

parms <- read.csv("../../Results/Figure2_Models_epi/Data_For_plots/value_functions/calibrated_parameters.csv")
parms$kappa_c <- 1
parms$kappa_l <- 1
parms$infolag <- 0 # How long are results delayed?
parms$test_quality <- 1 # How good are tests? When tests are perfect (t_q = 1), everyone perfectly knows their type. When tests are useless (t_q = 0), everyone has an equal chance of being told they are S, I, or R types. Accordingly useless tests result in the representative agents of each type playing the uniformly mixed strategies (1/3,1/3,1/3)%*%(S,I,R)
parms$eqm_concept <- "linear" # qre = solve for the quantal response equilibrium (takes a long time, only use on a big machine), linear = use linear interpolation
parms$compliance <- 1 # What proportion of people comply with policy? 

lag_sequence <- c(rep(14,length.out=60), rep(0,length.out=(final.time-60)))

infolag_list_benchmark <- list(lag_sequence = rep(0,length.out=final.time))
infolag_list <- list(lag_sequence = lag_sequence)

exog_parms <- parms

policy_list <- list(switch = 1, start_date = 35, end_date = 75, labor_supply_level = 0.65, state_dependent = 0)

#####
# Read in policy functions
#####

eqm_vpfn <- cbind(read.csv("../../Results/Figure2_Models_epi/Data_For_plots/value_functions/nextgen_eqm_vpfn.csv"),type="eqm")
opt_vpfn <- cbind(read.csv("../../Results/Figure2_Models_epi/Data_For_plots/value_functions/nextgen_planner_vpfn.csv"),type="plan")

# eqm_vpfn <- cbind(read.csv("../../Results/value_policy_functions/value_policy_functions__eqm_cheby_test.csv"),type="eqm")
# opt_vpfn <- cbind(read.csv("../../Results/value_policy_functions/value_policy_functions__planner_cheby_test.csv"),type="plan")

## Create grid and other necessary objects
S_gridlength <- length(unique(eqm_vpfn$S))
I_gridlength <- length(unique(eqm_vpfn$I))
R_gridlength <- length(unique(eqm_vpfn$R))
grid_pieces <- build_grid(S_gridlength, I_gridlength, R_gridlength)
grid_list <- grid_pieces[1:3]

grid_dfrm <- grid_pieces$dfrm
simplex_check <- rep(1,length.out=nrow(grid_dfrm))

setwd(paste0("../../Results/FigureX_compliance/", parms$eqm_concept, "_behavior"))

#### Block to read in pre-computed datasets

# dfrm_big_worst <- read_csv("Worst-case_dynamics.csv")
# dfrm_big_tqimp_8 <- read_csv("Improving test quality_dynamics.csv")
# dfrm_big_tqlow_85 <- read_csv("Improving test timeliness_dynamics.csv")
# dfrm_big_tqimp_85 <- read_csv("Improving test quality and timeliness_dynamics.csv")

#####
# Generate simulated time paths
#####

day_test_gets_perfect <- 75
max_test_quality <- 0.95

##### Perfect compliance

exog_parms$test_quality <- 1
exog_parms$infolag <- 0
exog_parms$compliance <- 1

eqm_fullcompliance <- generate_time_series(eqm_vpfn)
eqm_fullcompliance <- cbind(eqm_fullcompliance,infolag=0,type="eqm")
lockdown_fullcompliance <- generate_time_series(eqm_vpfn, policy_list=policy_list, fallback_solved_values=eqm_vpfn)
lockdown_fullcompliance <- cbind(lockdown_fullcompliance,infolag=0,type="lockdown")
opt_fullcompliance <- generate_time_series(opt_vpfn, fallback_solved_values=eqm_vpfn)
opt_fullcompliance <- cbind(opt_fullcompliance,infolag=0,type="plan")

dfrm_fullcompliance <- rbind(eqm_fullcompliance, lockdown_fullcompliance, opt_fullcompliance)
full_compliance_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_fullcompliance, "Full compliance")
full_compliance_figs_summary <- generate_recession_metrics_small_composite(dfrm_fullcompliance, "Full compliance")

full_compliance_figs_dynamics$infection + full_compliance_figs_dynamics$recession 
full_compliance_figs_summary$overall

##### Only 75% of people comply

exog_parms$compliance <- 0.75

eqm_partialcompliance <- generate_time_series(eqm_vpfn)
eqm_partialcompliance <- cbind(eqm_partialcompliance,infolag=0,type="eqm")
lockdown_partialcompliance <- generate_time_series(eqm_vpfn, policy_list=policy_list, fallback_solved_values=eqm_vpfn)
lockdown_partialcompliance <- cbind(lockdown_partialcompliance,infolag=0,type="lockdown")
opt_partialcompliance <- generate_time_series(opt_vpfn, fallback_solved_values=eqm_vpfn)
opt_partialcompliance <- cbind(opt_partialcompliance,infolag=0,type="plan")

dfrm_partialcompliance <- rbind(eqm_partialcompliance, lockdown_partialcompliance, opt_partialcompliance)
partial_compliance_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_partialcompliance, "Partial compliance, perfect information")
partial_compliance_figs_summary <- generate_recession_metrics_small_composite(dfrm_partialcompliance, "Partial compliance, perfect information")

partial_compliance_figs_summary$overall

##### Nobody complies
exog_parms$compliance <- 0

eqm_nocompliance <- generate_time_series(eqm_vpfn)
eqm_nocompliance <- cbind(eqm_nocompliance,infolag=0,type="eqm")

lockdown_nocompliance <- generate_time_series(eqm_vpfn, policy_list=policy_list, fallback_solved_values=eqm_vpfn)
lockdown_nocompliance <- cbind(lockdown_nocompliance,infolag=0,type="lockdown")

opt_nocompliance <- generate_time_series(opt_vpfn, fallback_solved_values=eqm_vpfn)
opt_nocompliance <- cbind(opt_nocompliance,infolag=0,type="plan")

dfrm_nocompliance <- rbind(eqm_nocompliance, lockdown_nocompliance, opt_nocompliance)
no_compliance_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_nocompliance, "Zero compliance, perfect information")
no_compliance_figs_summary <- generate_recession_metrics_small_composite(dfrm_nocompliance, "Zero compliance, perfect information")

no_compliance_figs_summary$cases + no_compliance_figs_summary$losses

##### Partial compliance with the "most realistic" information structure

quality_sequence <- c( seq(from=0.1, to=max_test_quality, length=day_test_gets_perfect), rep(1,length.out=(final.time-day_test_gets_perfect)))

lag_sequence <- c(rep(8,length.out=60), rep(5,length.out=(final.time-60)))
infolag_list_seq <- list(lag_sequence = lag_sequence, quality_sequence = quality_sequence)
exog_parms$compliance <- 1

eqm_partialcompliance2 <- generate_time_series(eqm_vpfn, infolag_list=infolag_list_seq)
eqm_partialcompliance2 <- cbind(eqm_partialcompliance2,infolag=0.85,type="eqm")
lockdown_partialcompliance2 <- generate_time_series(eqm_vpfn, policy_list=policy_list, fallback_solved_values=eqm_vpfn, infolag_list=infolag_list_seq)
lockdown_partialcompliance2 <- cbind(lockdown_partialcompliance2,infolag=0.85,type="lockdown")
opt_partialcompliance2 <- generate_time_series(opt_vpfn, fallback_solved_values=eqm_vpfn, infolag_list=infolag_list_seq)
opt_partialcompliance2 <- cbind(opt_partialcompliance2,infolag=0.85,type="plan")

dfrm_partialcompliance2 <- rbind(eqm_partialcompliance2, lockdown_partialcompliance2, opt_partialcompliance2)
partial_compliance2_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_partialcompliance2, "Partial compliance, improving information")
partial_compliance2_figs_summary <- generate_recession_metrics_small_composite(dfrm_partialcompliance2, "Partial compliance, improving information")


#####
# Scatterplot x-y
#####

metrics_nc_perf <- generate_just_metrics(dfrm_nocompliance)
metrics_pc_perf <- generate_just_metrics(dfrm_partialcompliance)
metrics_pc_2 <- generate_just_metrics(dfrm_partialcompliance2)

list_of_metrics <- list(metrics_nc_perf, metrics_pc_perf, metrics_pc_2)

loss_ratios <- c(maximum_gain, compare_ratios(list_of_metrics,"losses"))
case_ratios <- c(max_averted_cases, compare_ratios(list_of_metrics,"casesp100k"))
scenarios <- c("Baseline","Low compliance, perfect information","Partial compliance, perfect information","Partial compliance, improving information") 

ratio_comps_2 <- data.frame(scenarios = scenarios, loss_ratios = loss_ratios, case_ratios = case_ratios) %>% arrange(-loss_ratios) %>% mutate(perc_change_loss = loss_ratios/loss_ratios[1], perc_change_case = case_ratios/case_ratios[1])
ratio_comps_2

write.csv(ratio_comps_2, file=paste0("compliance_scenarios_summary_ratios.csv"))

friction_effect_summary <- ggplot(data = ratio_comps_2, aes(x=perc_change_loss, y=perc_change_case, group=scenarios, color=scenarios)) + 
	geom_point(size = 10) +
	theme_minimal() +
	# geom_text(aes(label=scenarios), position=position_dodge(height=0.9), vjust=-0.25) +
	xlim(c(0,1)) + ylim(c(0,1)) +
	geom_abline(intercept = 0, slope = 1, linetype="dashed") +
	labs(x="% of baseline economic benefits from targeting", y="% of baseline infection control benefits from targeting", title="Effects of compliance on policy effectiveness", color = "Scenarios") +
	scale_color_viridis(discrete=TRUE, option="turbo") +
            theme(strip.background = element_blank(),
                 strip.placement = "outside",
                 text = element_text(size=24)) 
friction_effect_summary


png(paste0("compliance_frictions_costs_",parms$eqm_concept,"_1.png"), width=1200, height=800)
friction_effect_summary
dev.off()

#####
# Draw figures
#####

## Perfect baseline
big_plot <- no_compliance_figs_summary$cases + no_compliance_figs_summary$losses + no_compliance_figs_dynamics$infection + no_compliance_figs_dynamics$recession + 
partial_compliance_figs_summary$cases + partial_compliance_figs_summary$losses + partial_compliance_figs_dynamics$infection + partial_compliance_figs_dynamics$recession + 
full_compliance_figs_summary$cases + full_compliance_figs_summary$losses + full_compliance_figs_dynamics$infection + full_compliance_figs_dynamics$recession + 
  plot_layout(ncol = 4, widths = c(0.75, 0.75, 1, 1), guides = "collect")
big_plot

png("../../Results/FigureX_compliance/compliance_fig_sketch__perfectbaseline.png", width=2000, height=1125)
big_plot
dev.off()


## Not perfect baseline
big_plot <- no_compliance_figs_summary$cases + no_compliance_figs_summary$losses + no_compliance_figs_dynamics$infection + no_compliance_figs_dynamics$recession + 
partial_compliance_figs_summary$cases + partial_compliance_figs_summary$losses + partial_compliance_figs_dynamics$infection + partial_compliance_figs_dynamics$recession + 
partial_compliance2_figs_summary$cases + partial_compliance2_figs_summary$losses + partial_compliance2_figs_dynamics$infection + partial_compliance2_figs_dynamics$recession  + 
  plot_layout(ncol = 4, widths = c(0.75, 0.75, 1, 1), guides = "collect")
big_plot

png("../../../Results/FigureX_compliance/compliance_fig_sketch__informationproblems.png", width=2000, height=1125)
big_plot
dev.off()



#####
# Assessing how low compliance has to go before targeted isolation produces identical costs to lockdown
#####

fc_lockdown_metrics <- generate_just_metrics(dfrm_fullcompliance) %>% filter(type=="lockdown") %>% rename(compliance = infolag)
compliance_rate_grid <- seq(0, 1, by=0.01)


lag_sequence <- c(rep(8,length.out=60), rep(5,length.out=(final.time-60)))
infolag_list_seq <- list(lag_sequence = lag_sequence, quality_sequence = quality_sequence)

varyingc_opt_metrics_list <- list()
varyingc_lockdown_metrics_list <- list()

for(compliance_rate_i in seq_along(compliance_rate_grid)) {
	compliance_grid_time <- proc.time()[3]
	exog_parms$compliance <- compliance_rate_grid[compliance_rate_i]

	opt_compliancegrid_temp <- generate_time_series(opt_vpfn, fallback_solved_values=eqm_vpfn, infolag_list=infolag_list_seq)
	opt_compliancegrid_temp <- cbind(opt_compliancegrid_temp,infolag=0,type="plan")

	lockdown_compliancegrid_temp <- generate_time_series(eqm_vpfn, policy_list=policy_list, fallback_solved_values=eqm_vpfn, infolag_list=infolag_list_seq)
	lockdown_compliancegrid_temp <- cbind(lockdown_compliancegrid_temp,infolag=0,type="lockdown")	

	opt_dfrm_temp <- rbind(eqm_fullcompliance, lockdown_fullcompliance, opt_compliancegrid_temp)
	varyingc_opt_metrics_list[[compliance_rate_i]] <- generate_just_metrics(opt_dfrm_temp) %>% filter(type=="coordinated") %>% rename(compliance = infolag) %>% mutate(compliance = exog_parms$compliance)

	lockdown_dfrm_temp <- rbind(eqm_fullcompliance, lockdown_compliancegrid_temp, opt_compliancegrid_temp)
	varyingc_lockdown_metrics_list[[compliance_rate_i]] <- generate_just_metrics(lockdown_dfrm_temp) %>% filter(type=="lockdown") %>% rename(compliance = infolag) %>% mutate(compliance = exog_parms$compliance)

	message("Completed compliance rate of ", compliance_rate_grid[compliance_rate_i],". Time taken: ", round( (proc.time()[3] - compliance_grid_time) ,3), " seconds.")
}

varyingc_lockdown_metrics_dfrm <- rbindlist(varyingc_lockdown_metrics_list)
varyingc_opt_metrics_dfrm <- rbindlist(varyingc_opt_metrics_list)
varyingc_metrics_dfrm <- rbind(varyingc_lockdown_metrics_dfrm, varyingc_opt_metrics_dfrm)

head(varyingc_metrics_dfrm)

write.csv(varyingc_metrics_dfrm,  file=paste0("breakeven_compliance_rate_assessment.csv"))
varyingc_metrics_dfrm <- read.csv(file=paste0("breakeven_compliance_rate_assessment.csv"))[,-1]

breakeven_compliance_rate_assessment_plot <- ggplot(data = varyingc_metrics_dfrm %>% filter(measure=="losses"), aes(x=compliance, group=type, linetype=type)) + geom_line(aes(y=values), size=1.5) + theme_bw()

breakeven_compliance_rate_assessment_plot 

png("breakeven_compliance_rate_assessment", width=1200, height=800)
breakeven_compliance_rate_assessment_plot
dev.off()

