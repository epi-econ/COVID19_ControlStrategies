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
# parms$compliance <- 1 # What proportion of people comply with policy? 

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

## Create grid and other necessary objects
S_gridlength <- length(unique(eqm_vpfn$S))
I_gridlength <- length(unique(eqm_vpfn$I))
R_gridlength <- length(unique(eqm_vpfn$R))
grid_pieces <- build_grid(S_gridlength, I_gridlength, R_gridlength)
grid_list <- grid_pieces[1:3]

grid_dfrm <- grid_pieces$dfrm
simplex_check <- rep(1,length.out=nrow(grid_dfrm))

setwd(paste0("../../Results/FigureX_infolag/", parms$eqm_concept, "_behavior"))

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

##### Perfect_baseline

exog_parms$test_quality <- 1
exog_parms$infolag <- 0

eqm_best <- generate_time_series(eqm_vpfn)
eqm_best <- cbind(eqm_best,type="eqm")
lockdown_best <- generate_time_series(eqm_vpfn, policy_list=policy_list)
lockdown_best <- cbind(lockdown_best,type="lockdown")
opt_best <- generate_time_series(opt_vpfn)
opt_best <- cbind(opt_best,type="plan")

# eqm_best_figs_dynamics <- generate_dynamics_figures_small(eqm_best, "Worst-case")
# lockdown_best_figs_dynamics <- generate_dynamics_figures_small(lockdown_best, "Worst-case")
# opt_best_figs_dynamics <- generate_dynamics_figures_small(opt_best, "Worst-case")

dfrm_big_best <- rbind(eqm_best, lockdown_best, opt_best)
best_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_best, "Limited and delayed testing")
best_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_best, "Limited and delayed testing")

# worst_figs_dynamics$infection
# worst_figs_dynamics$recession

best_figs_dynamics$infection | best_figs_dynamics$recession

##### Low test quality, 8 day test lag, compliant

exog_parms$test_quality <- 0.1
exog_parms$infolag <- 8

eqm_worst <- generate_time_series(eqm_vpfn)
eqm_worst <- cbind(eqm_worst,infolag=0.8,type="eqm")
lockdown_worst <- generate_time_series(eqm_vpfn, policy_list=policy_list)
lockdown_worst <- cbind(lockdown_worst,infolag=0.8,type="lockdown")
opt_worst <- generate_time_series(opt_vpfn)
opt_worst <- cbind(opt_worst,infolag=0.8,type="plan")

# eqm_worst_figs_dynamics <- generate_dynamics_figures_small(eqm_worst, "Worst-case")
# lockdown_worst_figs_dynamics <- generate_dynamics_figures_small(lockdown_worst, "Worst-case")
# opt_worst_figs_dynamics <- generate_dynamics_figures_small(opt_worst, "Worst-case")

dfrm_big_worst <- rbind(eqm_worst, lockdown_worst, opt_worst)
worst_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_worst, "Limited and delayed testing")
worst_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_worst, "Limited and delayed testing")

# worst_figs_dynamics$infection
# worst_figs_dynamics$recession

worst_figs_dynamics$infection | worst_figs_dynamics$recession

##### Improving test quality, 8 day test lag, compliant

quality_sequence <- c( seq(from=0.1, to=max_test_quality, length=day_test_gets_perfect), rep(1,length.out=(final.time-day_test_gets_perfect)))
infolag_list_seq <- list(lag_sequence = c(rep(8,length.out=final.time)), quality_sequence = quality_sequence)

eqm_tqimp_8 <- generate_time_series(eqm_vpfn, infolag_list=infolag_list_seq)
eqm_tqimp_8 <- cbind(eqm_tqimp_8,infolag=1.8,type="eqm")
lockdown_tqimp_8 <- generate_time_series(eqm_vpfn, infolag_list=infolag_list_seq, policy_list=policy_list)
lockdown_tqimp_8 <- cbind(lockdown_tqimp_8,infolag=1.8,type="lockdown")
opt_tqimp_8 <- generate_time_series(opt_vpfn, infolag_list=infolag_list_seq)
opt_tqimp_8 <- cbind(opt_tqimp_8,infolag=1.8,type="plan")

dfrm_big_tqimp_8 <- rbind(eqm_tqimp_8, lockdown_tqimp_8, opt_tqimp_8)
tqimp_8_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_tqimp_8, "Improving test quality")
tqimp_8_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_tqimp_8, "Improving test quality")

tqimp_8_figs_dynamics$infection | tqimp_8_figs_dynamics$recession


# tqimp_8_figs_dynamics$infection
# tqimp_8_figs_dynamics$recession

##### Low test quality, 8-5 improvement structure, compliant

quality_sequence <- c( rep(0.1,length.out=(final.time-day_test_gets_perfect)))
lag_sequence <- c(rep(8,length.out=60), rep(5,length.out=(final.time-60)))
infolag_list_seq <- list(lag_sequence = lag_sequence, quality_sequence = quality_sequence)

eqm_tqlow_85 <- generate_time_series(eqm_vpfn, infolag_list=infolag_list_seq)
eqm_tqlow_85 <- cbind(eqm_tqlow_85,infolag=0.85,type="eqm")
lockdown_tqlow_85 <- generate_time_series(eqm_vpfn, infolag_list=infolag_list_seq, policy_list=policy_list)
lockdown_tqlow_85 <- cbind(lockdown_tqlow_85,infolag=0.85,type="lockdown")
opt_tqlow_85 <- generate_time_series(opt_vpfn, infolag_list=infolag_list_seq)
opt_tqlow_85 <- cbind(opt_tqlow_85,infolag=0.85,type="plan")

dfrm_big_tqlow_85 <- rbind(eqm_tqlow_85, lockdown_tqlow_85, opt_tqlow_85)
tqlow_85_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_tqlow_85, "Improving test timeliness")
tqlow_85_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_tqlow_85, "Improving test timeliness")

# dfrm <- dfrm_big_tqlow_85

# tqlow_85_figs_dynamics$infection
# tqlow_85_figs_dynamics$recession

##### Improving test quality, 8-5 improvement structure, compliant

quality_sequence <- c( seq(from=0.1, to=max_test_quality, length=day_test_gets_perfect), rep(1,length.out=(final.time-day_test_gets_perfect)))
lag_sequence <- c(rep(8,length.out=60), rep(5,length.out=(final.time-60)))
infolag_list_seq <- list(lag_sequence = lag_sequence, quality_sequence = quality_sequence)

eqm_tqimp_85 <- generate_time_series(eqm_vpfn, infolag_list=infolag_list_seq)
eqm_tqimp_85 <- cbind(eqm_tqimp_85,infolag=1.85,type="eqm")
lockdown_tqimp_85 <- generate_time_series(eqm_vpfn, infolag_list=infolag_list_seq, policy_list=policy_list)
lockdown_tqimp_85 <- cbind(lockdown_tqimp_85,infolag=1.85,type="lockdown")
opt_tqimp_85 <- generate_time_series(opt_vpfn, infolag_list=infolag_list_seq)
opt_tqimp_85 <- cbind(opt_tqimp_85,infolag=1.85,type="plan")

dfrm_big_tqimp_85 <- rbind(eqm_tqimp_85, lockdown_tqimp_85, opt_tqimp_85)
tqimp_85_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_tqimp_85, "Improving test quality and delays")
tqimp_85_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_tqimp_85, "Improving test quality and delays")

# tqimp_85_figs_dynamics$infection
# tqimp_85_figs_dynamics$recession

#####
# Time path figures
#####

big_plot <- worst_figs_summary$overall + worst_figs_dynamics$infection + worst_figs_dynamics$recession + 
tqimp_8_figs_summary$overall + tqimp_8_figs_dynamics$infection + tqimp_8_figs_dynamics$recession + 
tqlow_85_figs_summary$overall + tqlow_85_figs_dynamics$infection + tqlow_85_figs_dynamics$recession + 
tqimp_85_figs_summary$overall + tqimp_85_figs_dynamics$infection + tqimp_85_figs_dynamics$recession + 
  plot_layout(ncol = 3, widths = c(1.5, 1, 1), guides = "collect")

big_plot <- worst_figs_summary$cases + worst_figs_summary$losses + worst_figs_dynamics$infection + worst_figs_dynamics$recession + 
tqlow_85_figs_summary$cases + tqlow_85_figs_summary$losses + tqlow_85_figs_dynamics$infection + tqlow_85_figs_dynamics$recession + 
tqimp_8_figs_summary$cases + tqimp_8_figs_summary$losses + tqimp_8_figs_dynamics$infection + tqimp_8_figs_dynamics$recession + 
tqimp_85_figs_summary$cases + tqimp_85_figs_summary$losses + tqimp_85_figs_dynamics$infection + tqimp_85_figs_dynamics$recession + 
  plot_layout(ncol = 4, widths = c(0.75, 0.75, 1, 1), guides = "collect")
# big_plot


#####
# Write out time path figures
#####

worst_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_worst, "Worst-case")
worst_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_worst, "Worst-case")

tqimp_8_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_tqimp_8, "Improving test quality")
tqimp_8_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_tqimp_8, "Improving test quality")

tqlow_85_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_tqlow_85, "Improving test timeliness")
tqlow_85_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_tqlow_85, "Improving test timeliness")

tqimp_85_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_tqimp_85, "Improving test quality and timeliness")
tqimp_85_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_tqimp_85, "Improving test quality and timeliness")

big_plot <- worst_figs_summary$cases + worst_figs_summary$losses + worst_figs_dynamics$infection + worst_figs_dynamics$recession + 
tqimp_8_figs_summary$cases + tqimp_8_figs_summary$losses + tqimp_8_figs_dynamics$infection + tqimp_8_figs_dynamics$recession + 
tqimp_85_figs_summary$cases + tqimp_85_figs_summary$losses + tqimp_85_figs_dynamics$infection + tqimp_85_figs_dynamics$recession + 
  plot_layout(ncol = 4, widths = c(0.75, 0.75, 1, 1), guides = "collect")
big_plot

png(paste0("infolag_fig_sketch_",parms$eqm_concept,"_cheby.png"), width=2000, height=1125)
big_plot
dev.off()

#####
# Scatterplot x-y
#####

metrics_worst <- generate_just_metrics(dfrm_big_worst)
metrics_tqimp_8 <- generate_just_metrics(dfrm_big_tqimp_8)
metrics_tqlow_85 <- generate_just_metrics(dfrm_big_tqlow_85)
metrics_tqimp_85 <- generate_just_metrics(dfrm_big_tqimp_85)

list_of_metrics <- list(metrics_worst, metrics_tqimp_8, metrics_tqimp_85)

loss_ratios <- c(maximum_gain, compare_ratios(list_of_metrics,"losses"))
case_ratios <- c(max_averted_cases, compare_ratios(list_of_metrics,"casesp100k"))
scenarios <- c("Frictionless baseline","Limited and delayed testing","Improving test quality","Improving test quality and delays") 

ratio_comps_2 <- data.frame(scenarios = scenarios, loss_ratios = loss_ratios, case_ratios = case_ratios) %>% arrange(-loss_ratios) %>% mutate(perc_change_loss = loss_ratios/loss_ratios[1], perc_change_case = case_ratios/case_ratios[1])
ratio_comps_2

write.csv(ratio_comps_2, file=paste0("infolag_scenarios_summary_ratios.csv"))


friction_effect_summary_2 <- ggplot(data = ratio_comps_2, aes(x=perc_change_loss, y=perc_change_case, group=scenarios, color=scenarios)) + 
	geom_point(size = 10) +
	theme_minimal() +
	# geom_text(aes(label=scenarios), position=position_dodge(height=0.9), vjust=-0.25) +
	xlim(c(0,1)) + ylim(c(0,1)) +
	geom_abline(intercept = 0, slope = 1, linetype="dashed") +
	labs(x="% of baseline economic benefits from targeting", y="% of baseline infection control benefits from targeting", title="Effects of frictions on policy effectiveness", color = "Scenarios") +
	scale_color_viridis(discrete=TRUE, option="turbo") +
            theme(strip.background = element_blank(),
                 strip.placement = "outside",
                 text = element_text(size=24)) 
friction_effect_summary_2


# png(paste0("info_frictions_costs_",parms$eqm_concept,"_1.png"), width=1200, height=800)
# friction_effect_summary_1
# dev.off()

png(paste0("info_frictions_costs_",parms$eqm_concept,"_2.png"), width=1200, height=800)
friction_effect_summary_2
dev.off()

#################################################################################################
# Diagnosing the mixing probabilities
#################################################################################################

# mixing_probs_diagnosis <- function(dataset) {

# 	mixing_probs_data <- dataset %>% select(starts_with("l_")|"type"|"time")
# 	head(mixing_probs_data)

# 	colnames(mixing_probs_data) <- NULL

# 	mixing_probs_data_long <- rbind( 
# 			(cbind(mixing_probs_data[,c(1:3,10:11)],agent_type="S")),
# 			(cbind(mixing_probs_data[,c(4:6,10:11)],agent_type="I")),
# 			(cbind(mixing_probs_data[,c(7:9,10:11)],agent_type="R")))
# 	colnames(mixing_probs_data_long) <- c("S_prob","I_prob","R_prob","policy_type","time","agent_type")

# 	mixing_probs_data_longer <- pivot_longer(mixing_probs_data_long, contains("_prob"), values_to="probability", names_to="choice")
# 	head(mixing_probs_data_longer)

# 	probs_plotlist <- list()
# 	ls_plotlist <- list()
# 	for(i in seq_along(unique(mixing_probs_data_longer$policy_type))){
# 		selection <- unique(mixing_probs_data_longer$policy_type)[i]
		
# 		probs_plotlist[[i]] <- ggplot(data = (mixing_probs_data_longer %>% filter(policy_type==selection)), aes(x=time, y=probability, group=choice, color=choice)) +
# 			geom_line(size=1) +
# 			geom_vline(xintercept = day_test_gets_perfect, linetype="dashed") +
# 			facet_wrap(~agent_type, nrow=1,labeller = as_labeller(c(S = "Type S choices", I = "Type I choices", R = "Type R choices"))) +
# 			labs(x="Day", y="QRE choice probability", title=paste0("Policy type: ",selection)) +
# 			theme_bw()

# 		labor_supplies_data <- dataset %>% select(c(time,S,I,R,labor_S,labor_I,labor_R,type)) %>% 
# 		mutate(S_prob = S*labor_S, I_prob = I*labor_I, R_prob = R*labor_R) %>%
# 		select("time"|contains("_prob"),"type") %>%
# 		pivot_longer(contains("_prob"), values_to="labor_supply", names_to="hours") %>%
# 		filter(type==selection)

# 		ls_plotlist[[i]] <- ggplot(data = labor_supplies_data, aes(x=time, y=labor_supply, group=hours, color=hours)) + 
# 			geom_line(size=1) +
# 			theme_bw() + ggtitle(paste0("Labor supplies: ",selection)) +
# 			labs(x="Day", y="Hours worked", color="choice") +
# 			geom_vline(xintercept = day_test_gets_perfect, linetype="dashed")
# 	}

# 	mixing_plot <- (probs_plotlist[[1]] | ls_plotlist[[1]]) / (probs_plotlist[[2]] | ls_plotlist[[2]]) / (probs_plotlist[[3]] | ls_plotlist[[3]])  + 
# 	  plot_layout(nrow = 3, guides = "collect")

# 	return(mixing_plot)

# }

# tqimp_85_md <- mixing_probs_diagnosis(dfrm_big_tqimp_85)