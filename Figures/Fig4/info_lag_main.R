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

# setwd("~/Dropbox/Corona_epicon/Code/methods_paper/")

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
parms$caseload_wt <- 0
parms$caseload_type <- "linear"
parms$infolag <- 0 # How long are results delayed?
parms$test_quality <- 1 # How good are tests? When tests are perfect (t_q = 1), everyone perfectly knows their type. When tests are useless (t_q = 0), everyone has an equal chance of being told they are S, I, or R types. Accordingly useless tests result in the representative agents of each type playing the uniformly mixed strategies (1/3,1/3,1/3)%*%(S,I,R)
parms$eqm_concept <- "linear" # qre = solve for the quantal response equilibrium (takes a long time, only use on a big machine), linear = use linear interpolation
# parms$compliance <- 1 # What proportion of people comply with policy? 
parms$final.time <- final.time

lag_sequence <- c(rep(14,length.out=60), rep(0,length.out=(final.time-60)))

infolag_list_benchmark <- list(lag_sequence = rep(0,length.out=final.time))
infolag_list <- list(lag_sequence = lag_sequence)

exog_parms <- parms

policy_list <- list(switch = 1, start_date = 35, end_date = 75, labor_supply_level = 0.65, state_dependent = 0)

label <- "8_5_3_lag"

#####
# Read in policy functions
#####

eqm_vpfn <- cbind(read.csv("../../Results/Figure2_Models_epi/Data_For_plots/value_functions/nextgen_eqm_vpfn.csv"),type="eqm")
opt_vpfn <- cbind(read.csv("../../Results/Figure2_Models_epi/Data_For_plots/value_functions/nextgen_planner_vpfn.csv"),type="plan")

# eqm_vpfn <- cbind(read.csv("../../Results/FigureX_infolag/qre_behavior/vpfn_eqm_cheby_16.csv"),type="eqm")
# opt_vpfn <- cbind(read.csv("../../Results/FigureX_infolag/qre_behavior/vpfn_plan_cheby_16.csv"),type="plan")

# eqm_vpfn <- cbind(read.csv("../../Results/value_policy_functions/value_policy_functions__eqm_cheby_test_24.csv"),type="eqm")
# opt_vpfn <- cbind(read.csv("../../Results/value_policy_functions/value_policy_functions__planner_cheby_test_24.csv"),type="plan")

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
max_test_quality <- 0.75

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
dfrm_big_best$infolag <- exog_parms$infolag
best_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_best, "Limited and delayed testing", end_time=final.time)
best_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_best, "Limited and delayed testing", end_time=final.time)

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
dfrm_big_best$infolag <- exog_parms$infolag
worst_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_worst, "Limited and delayed testing", end_time=final.time)
worst_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_worst, "Limited and delayed testing", end_time=final.time)

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
dfrm_big_best$infolag <- exog_parms$infolag
tqimp_8_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_tqimp_8, "Improving test quality", end_time=final.time)
tqimp_8_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_tqimp_8, "Improving test quality", end_time=final.time)

tqimp_8_figs_dynamics$infection | tqimp_8_figs_dynamics$recession


# tqimp_8_figs_dynamics$infection
# tqimp_8_figs_dynamics$recession

##### Low test quality, 8-5 (8-5-3) improvement structure, compliant

# quality_sequence <- c( rep(0.1,length.out=(final.time-day_test_gets_perfect)))
quality_sequence <- c( rep(0.1,length.out=final.time))
# lag_sequence <- c(rep(8,length.out=60), rep(5,length.out=(final.time-60))) # 8-5 lag structure
lag_sequence <- c(
		rep(8,length.out=60), 
		rep(5,length.out=(75-60)),
		rep(3, length.out=(final.time-75))
		 )
infolag_list_seq <- list(lag_sequence = lag_sequence, quality_sequence = quality_sequence)

eqm_tqlow_85 <- generate_time_series(eqm_vpfn, infolag_list=infolag_list_seq)
eqm_tqlow_85 <- cbind(eqm_tqlow_85,infolag=0.85,type="eqm")
lockdown_tqlow_85 <- generate_time_series(eqm_vpfn, infolag_list=infolag_list_seq, policy_list=policy_list)
lockdown_tqlow_85 <- cbind(lockdown_tqlow_85,infolag=0.85,type="lockdown")
opt_tqlow_85 <- generate_time_series(opt_vpfn, infolag_list=infolag_list_seq)
opt_tqlow_85 <- cbind(opt_tqlow_85,infolag=0.85,type="plan")

dfrm_big_tqlow_85 <- rbind(eqm_tqlow_85, lockdown_tqlow_85, opt_tqlow_85)
dfrm_big_best$infolag <- exog_parms$infolag
tqlow_85_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_tqlow_85, "Improving test timeliness", end_time=final.time)
tqlow_85_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_tqlow_85, "Improving test timeliness", end_time=final.time)

# dfrm <- dfrm_big_tqlow_85

# tqlow_85_figs_dynamics$infection
# tqlow_85_figs_dynamics$recession

##### Improving test quality, 8-5 improvement structure, compliant

quality_sequence <- c( seq(from=0.1, to=max_test_quality, length=day_test_gets_perfect), rep(1,length.out=(final.time-day_test_gets_perfect)))
# lag_sequence <- c(rep(8,length.out=60), rep(5,length.out=(final.time-60)))
lag_sequence <- c(
		rep(8,length.out=60), 
		rep(5,length.out=(75-60)),
		rep(3, length.out=(final.time-75))
		 )
infolag_list_seq <- list(lag_sequence = lag_sequence, quality_sequence = quality_sequence)

eqm_tqimp_85 <- generate_time_series(eqm_vpfn, infolag_list=infolag_list_seq)
eqm_tqimp_85 <- cbind(eqm_tqimp_85,infolag=1.85,type="eqm")
lockdown_tqimp_85 <- generate_time_series(eqm_vpfn, infolag_list=infolag_list_seq, policy_list=policy_list)
lockdown_tqimp_85 <- cbind(lockdown_tqimp_85,infolag=1.85,type="lockdown")
opt_tqimp_85 <- generate_time_series(opt_vpfn, infolag_list=infolag_list_seq)
opt_tqimp_85 <- cbind(opt_tqimp_85,infolag=1.85,type="plan")

dfrm_big_tqimp_85 <- rbind(eqm_tqimp_85, lockdown_tqimp_85, opt_tqimp_85)
dfrm_big_best$infolag <- exog_parms$infolag
tqimp_85_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_tqimp_85, "Improving test quality and delays", end_time=final.time)
tqimp_85_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_tqimp_85, "Improving test quality and delays", end_time=final.time)

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

worst_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_worst, "Worst-case", end_time=final.time)
worst_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_worst, "Worst-case", end_time=final.time)

tqimp_8_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_tqimp_8, "Improving test quality", end_time=final.time)
tqimp_8_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_tqimp_8, "Improving test quality", end_time=final.time)

tqlow_85_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_tqlow_85, "Improving test timeliness", end_time=final.time)
tqlow_85_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_tqlow_85, "Improving test timeliness", end_time=final.time)

tqimp_85_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_tqimp_85, "Improving test quality and timeliness", end_time=final.time)
tqimp_85_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_tqimp_85, "Improving test quality and timeliness", end_time=final.time)

big_plot <- worst_figs_summary$cases + worst_figs_summary$losses + worst_figs_dynamics$infection + worst_figs_dynamics$recession + 
tqimp_8_figs_summary$cases + tqimp_8_figs_summary$losses + tqimp_8_figs_dynamics$infection + tqimp_8_figs_dynamics$recession + 
tqimp_85_figs_summary$cases + tqimp_85_figs_summary$losses + tqimp_85_figs_dynamics$infection + tqimp_85_figs_dynamics$recession + 
  plot_layout(ncol = 4, widths = c(0.75, 0.75, 1, 1), guides = "collect")
big_plot

png(paste0("infolag_fig_sketch_",parms$eqm_concept,"_",label,".png"), width=2000, height=1125)
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

write.csv(ratio_comps_2, file=paste0("infolag_scenarios_summary_ratios_",label,".csv"))


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

png(paste0("info_frictions_costs_",parms$eqm_concept,"_",label,".png"), width=1200, height=800)
friction_effect_summary_2
dev.off()

#####
# Test out a whole mess of different delays
#####

list_mess_data <- list()
list_mess_figs <- list()

day_test_gets_perfect <- 1
max_test_quality <- 1

##### A mess

for(quality in seq(0,1,length.out=25)){
	exog_parms$test_quality <- quality

	for(delay in 1:15) {
		exog_parms$infolag <- delay

		eqm_mess <- generate_time_series(eqm_vpfn)
		eqm_mess <- cbind(eqm_mess,type=paste0("eqm_lag_",delay,"_quality_",exog_parms$test_quality))
		lockdown_mess <- generate_time_series(eqm_vpfn, policy_list=policy_list)
		lockdown_mess <- cbind(lockdown_mess,type=paste0("lockdown_lag_",delay,"_quality_",exog_parms$test_quality))
		opt_mess <- generate_time_series(opt_vpfn)
		opt_mess <- cbind(opt_mess,type=paste0("plan_lag_",delay,"_quality_",exog_parms$test_quality))

		dfrm_big_mess <- rbind(eqm_mess, lockdown_mess, opt_mess)
		dfrm_big_mess$infolag <- exog_parms$infolag
		# mess_figs_dynamics <- generate_dynamics_figures_small_composite(dfrm_big_mess, "Limited and delayed testing")
		# mess_figs_summary <- generate_recession_metrics_small_composite(dfrm_big_mess, "Limited and delayed testing")

		# figs_list <- list(dynamics=mess_figs_dynamics, summary=mess_figs_summary)

		list_mess_data[[delay]] <- dfrm_big_mess
		# list_mess_figs[[delay]] <- figs_list


		# best_figs_dynamics$infection | best_figs_dynamics$recession

	}
	message(paste0("Quality ", round(quality,3), " complete"))
}

dfrm_mess <- rbindlist(list_mess_data)

write.csv(dfrm_mess, "delay_sims.csv")

dfrm_mess <- read.csv("delay_sims.csv")[,-1]

dfrm_mess$type <- word(dfrm_mess$type, start=1L, sep=fixed("_"))

dfrm_mess_sum <- dfrm_mess %>% 
				mutate(weighted_labor_S = S*labor_S) %>%
				mutate(weighted_labor_I = I*labor_I) %>%
				mutate(weighted_labor_R = R*labor_R) %>% 
				mutate(total_contacts = consumption_contacts + labor_contacts + other_contacts) %>%
				mutate(type = recode(type, eqm = "decentralized", lockdown = "lockdown", plan = "coordinated")) %>%
				mutate(SI_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_I, l_S = labor_S, l_I = labor_I) )) %>%
				mutate(SR_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_R, l_S = labor_S, l_I = labor_R) )) %>%
				mutate(RI_contacts = contacts(parms, choices = data.frame(c_S = consumption_R, c_I = consumption_I, l_S = labor_R, l_I = labor_I) )) %>%
				mutate(SS_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_S, l_S = labor_S, l_I = labor_S) )) %>%
				mutate(RR_contacts = contacts(parms, choices = data.frame(c_S = consumption_R, c_I = consumption_R, l_S = labor_R, l_I = labor_R) )) %>%
				mutate(prob_contact_I = (SI_contacts + RI_contacts)/(SI_contacts + RI_contacts + SR_contacts + SS_contacts + RR_contacts)) %>%
				mutate(prob_contact_I_weighted = (SI_contacts*S*I + RI_contacts*R*I)/(SI_contacts*S*I + RI_contacts*R*I + SR_contacts*S*R + SS_contacts*S*S + RR_contacts*R*R)) %>%
				mutate(daily_new_cases_from_consumption = parms$tau*parms$rho_c*consumption_S*consumption_I*S*I) %>%
				mutate(daily_new_cases_from_labor = parms$tau*parms$rho_l*labor_S*labor_I*S*I) %>%
				mutate(daily_new_unavoidable_cases = parms$tau*parms$rho_o*S*I) #%>%
				#mutate(quality = str_split(type, "_", simplify=TRUE)[,5]) # re-run this across multiple qualities

dfrm_mess_wide <- pivot_wider(dfrm_mess_sum[,c("time","type","aggregate_consumption","prob_contact_I_weighted","I","R","D","infolag")], id_cols=c("time","type","infolag"), names_from=type, values_from=c(aggregate_consumption, prob_contact_I_weighted, I, R, D), names_sep="_")


  for(i in seq_along(unique(dfrm_mess_sum[,"type"]))) {
      core_name <- unique(dfrm_mess_sum[,"type"])[i]
      new_colname <- paste0("PV_",core_name,"_losses")
      cons_colname <- paste0("aggregate_consumption_",core_name)
      print(cons_colname)
      dfrm_mess_wide <- dfrm_mess_wide %>% mutate(!!new_colname := (58000/365 - !!rlang::sym(cons_colname))*discount_factor^time)
  }


temp <- list() 

for(i in seq_along(unique(dfrm_mess_sum[,"type"]))) {
    core_name <- unique(dfrm_mess_sum[,"type"])[i]
    loss_colname <- paste0(core_name,"_losses")
    case_colname <- paste0(core_name,"_casesp100k")
    peak_colname <- paste0(core_name,"_peakcases")
    vol_colname <- paste0(core_name,"_volatility")

    loss_colname_input <- paste0("PV_",core_name,"_losses")

    I_name <- paste0("I_",core_name)
    R_name <- paste0("R_",core_name)
    D_name <- paste0("D_",core_name)

    temp_wide <- dfrm_mess_wide %>% group_by(infolag) %>% 
        summarise(!!loss_colname := round(sum( !!rlang::sym(loss_colname_input), na.rm=TRUE )),
                  !!case_colname := round(max(!!rlang::sym(D_name) + !!rlang::sym(I_name) + !!rlang::sym(R_name))*100000 ,0 ),
                  !!peak_colname := max(!!rlang::sym(I_name)),
                  !!vol_colname := (sd((!!rlang::sym(I_name)))))

    temp[[i]] <- pivot_longer(temp_wide, cols=c(
        !!rlang::sym(loss_colname), !!rlang::sym(case_colname), !!rlang::sym(peak_colname), !!rlang::sym(vol_colname), 
        ), names_to = c("type","measure"), names_sep="_", values_to = "values")
}

bar_data <- rbindlist(temp) %>% 
				mutate(type = recode(type, decentralized = "Voluntary\nisolation", lockdown = "Blanket\nlockdown", coordinated = "Targeted\nisolation"))

losses <- ggplot((bar_data %>% filter(measure=="losses")), aes(fill=type, x=infolag, y=values)) + 
	geom_bar(position="dodge", stat="identity") +
	labs(x="Information delay (days)", y="Economic loss ($/person)", fill="Policy type") +
	theme_bw() +
            theme(strip.background = element_blank(),
                 strip.placement = "outside",
                 text = element_text(size=24))

cases <- ggplot((bar_data %>% filter(measure=="casesp100k")), aes(fill=type, x=infolag, y=values)) + 
	geom_bar(position="dodge", stat="identity") +
	labs(x="Information delay (days)", y="Total cases (cases/100k)", fill="Policy type") +
	theme_bw() +
            theme(strip.background = element_blank(),
                 strip.placement = "outside",
                 text = element_text(size=24))

peaks <- ggplot((bar_data %>% filter(measure=="peakcases")), aes(fill=type, x=infolag, y=values)) + 
	geom_bar(position="dodge", stat="identity") +
	labs(x="Information delay (days)", y="Peak cases (% of population)", fill="Policy type") +
	theme_bw() +
            theme(strip.background = element_blank(),
                 strip.placement = "outside",
                 text = element_text(size=24))

vols <- ggplot((bar_data %>% filter(measure=="volatility")), aes(fill=type, x=infolag, y=values)) + 
	geom_bar(position="dodge", stat="identity") +
	labs(x="Information delay (days)", y="Daily caseload sd (% of population)", fill="Policy type") +
	theme_bw() +
            theme(strip.background = element_blank(),
                 strip.placement = "outside",
                 text = element_text(size=24))

combined_plot <- losses + cases + peaks + vols + plot_layout(guides='collect', ncol = 2)

combined_plot


png("info_delay_summary.png", width=2000, height=1125)
combined_plot
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