##############################################################################
# Simulation of epicon model
# Script layout:

##############################################################################

##### Load packages and scripts, set options, accept command line arguments

PD_range <- seq(from=0.005,to=0.1,length.out=10)

for(sensitivity in seq_along(PD_range)) {

library(ggplot2)
library(cowplot)
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

options(width=100)
options(scipen=100)

enableJIT(3) 

source("functions.R")
source("epicon_calibration_functions.R")

scriptargs = commandArgs(trailingOnly = TRUE)

#####
# A: Create grid, set non-utility structural parameters prior to calibration
#####

ncores = 64

## initial conditions: placed here so that they can be used to generate the grid
total_population = 331002651				# 331002651 is the US population
I_0 = 331.002651/total_population			# set to get 1-in-a-million as an initial condition
S_0 = 1 - I_0
R_0 = 0
SIR_init = data.frame(S=S_0, I=I_0, R=R_0)

## interactive run arguments
if(length(scriptargs)==0) {
	S_gridlength = 15
	I_gridlength = 15
	R_gridlength = 15
	grid_pieces = build_grid(S_gridlength, I_gridlength, R_gridlength)
	hot_start = 0

	social = 0
	problemtype = "decentralized"

	# structural infection/contact parameters
	daily_cons_contacts = 4.2877 		# avg daily contacts at consumption activities
	daily_work_contacts = 3.1625		# avg daily contacts at labor activities
	daily_other_contacts = 5			# avg daily unavoidable contacts
	R0_target_parameter = 2.6
	duration = 5.1						# duration of infectious period

	## economic calibration targets
	Cbm = (58000/365)	# daily annuitized consumption flow in dollars
	risk_aversion = 0.1

	calibrate = 0

	scenario_label = "interactive_run"

	phi = 0.8554632		# labor productivity of infected type, others are 1. The calculation is shown in the SI Appendix.

	precision = 0

	# target_R0_or_tau = "R0"
	# R0_target_parameter = 2.6
	target_R0_or_tau = "tau"
	R0_target_parameter = 0.04094745

	kappa_c = 1
	kappa_l = 1
}
## command line arguments
if(length(scriptargs)>0) {
	print(scriptargs)
	ncores = as.numeric(scriptargs[1])
	S_gridlength = as.numeric(scriptargs[2])
	I_gridlength = as.numeric(scriptargs[2])
	R_gridlength = as.numeric(scriptargs[2])
	grid_pieces = build_grid(S_gridlength, I_gridlength, R_gridlength)
	hot_start = as.numeric(scriptargs[3])

	social = ifelse(as.character(scriptargs[4])=="social", social <- 1, social <- 0)
	ifelse(as.character(scriptargs[5])=="planner", problemtype <- "planner", problemtype <- "decentralized")

	daily_cons_contacts = as.numeric(scriptargs[6])		# avg daily contacts at consumption activities
	daily_work_contacts = as.numeric(scriptargs[7])		# avg daily contacts at labor activities
	daily_other_contacts = as.numeric(scriptargs[8])	# avg daily unavoidable contacts
	R0_target_parameter = as.numeric(scriptargs[9])		# determines value of tau_parm if targeting R0, or the value of R0 if targeting tau (see target_R0_or_tau)
	duration = as.numeric(scriptargs[10])		# duration of infectious period

	Cbm = as.numeric(scriptargs[11])/365	# daily annuitized consumption flow in dollars
	risk_aversion = as.numeric(scriptargs[12])		# risk aversion; U_D is calibrated to adjust to this

	calibrate = as.numeric(scriptargs[13])	# binary flag for recalibrating utility parameters or not. set to 1 if this is the first time running the script with these parameters.

	scenario_label = paste0(as.character(scriptargs[14]),"_PD_",round(PD_range[sensitivity],4))	# scenario label

	target_R0_or_tau = as.character(scriptargs[15])		# should we target R0 or tau? by default will assume we're targeting R0

	precision = as.numeric(scriptargs[16])		# level of precision in VFI

	phi = as.numeric(scriptargs[17])		# labor productivity of infected type, others are 1. The calculation is shown in the SI Appendix.
	if(length(scriptargs)>17) {
		kappa_c = as.numeric(scriptargs[18])
		kappa_l = as.numeric(scriptargs[19])
	}
	if(length(scriptargs)<=17) {
		kappa_c = 1
		kappa_l = 1
	}
}

Lbm = 0.3333*24		# hours of day spent working
Lbar = 12 			# hours of day available for leisure
plannertype = "weighted" 
measurezero = 1
wage = Cbm/Lbm 		# hourly wage
policy = 0
nonlabor_income = 0
VSL = 1e07			# Value of a statistical life. Set to $10m.

discount_rate = 0.04								# annual discount rate
discount_factor = (1/(1+discount_rate))^(1/365)		# daily discount factor

final.time = 548	# number of periods for time series projection

grid_list = grid_pieces[1:3]
grid_dfrm = grid_pieces$dfrm

simplex_check = rep(1,length.out=nrow(grid_dfrm))

rho_c = daily_cons_contacts/(Cbm^(2*kappa_c))	# daily_cons_contacts = rho_c*Cbm*Cbm
rho_l = daily_work_contacts/(Lbm^(2*kappa_l))	# daily_work_contacts = rho_l*Lbm*Lbm
rho_o = daily_other_contacts		# O*O := 1 (normalization)

est_cfr = PD_range[sensitivity]						# conditional on leaving I, 5% chance end up in D; Pr(I->D)
proportion_removed = 1/duration		# 1 - (1 - (1/duration))^7
pi_r = (1-est_cfr)*proportion_removed
pi_d = est_cfr*proportion_removed

if(target_R0_or_tau=="R0"){
	target_R0 = R0_target_parameter
	tau_parm = target_R0/(duration*(daily_cons_contacts+daily_work_contacts+daily_other_contacts))
}
if(target_R0_or_tau=="tau") {
	tau_parm = R0_target_parameter
	target_R0 = tau_parm*(duration*(daily_cons_contacts+daily_work_contacts+daily_other_contacts))
}
#####
# B. Calibrate utility parameters
#####
source("calibrate_utility_parameters.R", print.eval=TRUE)

s = best_parms$s 				# elasticity of substitution between leisure and consumption
alpha_U = best_parms$alpha_U	# CES share parameter
gamma_c = best_parms$gamma_c 	# consumption contact utility parameter
gamma_l = best_parms$gamma_l 	# labor contact utility parameter
U_D = best_parms$U_D

l_S_guess = values_at_best_parms$labor_supply
l_I_guess = values_at_best_parms$labor_supply_I
l_R_guess = values_at_best_parms$labor_supply

#####
# C. Generate value function guess
#####

## create clusters
cl = makeCluster(ncores)
registerDoParallel(cl)

exog_parms = data.frame(gamma_c=gamma_c, gamma_l=gamma_l, rho_c=rho_c, rho_l=rho_l, rho_o=rho_o, pi_r=pi_r, pi_d=pi_d, discount_factor=discount_factor, phi=phi, w=wage, final.time=final.time, s=s, risk_aversion=risk_aversion, gamma = duration, tau_parm=tau_parm, R0=target_R0, alpha_U=alpha_U, precision=precision, udeath = U_D, Cbm = Cbm, Lbm = Lbm, Lbar = Lbar, problemtype = problemtype, plannertype = plannertype, kappa_c = kappa_c, kappa_l = kappa_l)

fwrite(exog_parms, file = paste0("../../Results/value_policy_functions/calibrated_parameters_",scenario_label,".csv"))
fwrite(exog_parms, file = paste0("../../Results/value_policy_functions/exog_parms.csv"))
# fwrite(exog_parms, file = paste0("../../Results/value_policy_functions/exog_parms_main.csv"))

print(grid_pieces[1])
print(grid_pieces[2])
print(grid_pieces[3])
print(sum(simplex_check))

print(scenario_label)

message("Using ", ncores, " cores.")

print(exog_parms)

print(U_D)

## continuation value guess
if(hot_start == 0) {
	ss_args_s = data.frame(c=wage*l_S_guess, l=l_S_guess, C=wage*l_S_guess, L=l_S_guess, gamma_c=gamma_c, gamma_l=gamma_l, s=s, risk_aversion=risk_aversion, phi=phi, rho_c=rho_c, rho_l=rho_l, rho_o=rho_o, pi_r=pi_r, pi_d=pi_d, alpha_U=alpha_U, Lbar=Lbar, kappa_c = kappa_c, kappa_l = kappa_l)

	ss_args_i = data.frame(c=phi*wage*l_I_guess, l=l_I_guess, C=wage*l_S_guess, L=l_S_guess, gamma_c=gamma_c, gamma_l=gamma_l, s=s, risk_aversion=risk_aversion, phi=phi, rho_c=rho_c, rho_l=rho_l, rho_o=rho_o, pi_r=pi_r, pi_d=pi_d, alpha_U=alpha_U, Lbar=Lbar, kappa_c = kappa_c, kappa_l = kappa_l)

	ss_args_r = data.frame(c=wage*l_R_guess, l=l_R_guess, C=wage*l_R_guess, L=l_R_guess, gamma_c=gamma_c, gamma_l=gamma_l, s=s, risk_aversion=risk_aversion, phi=phi, rho_c=rho_c, rho_l=rho_l, rho_o=rho_o, pi_r=pi_r, pi_d=pi_d, alpha_U=alpha_U, Lbar=Lbar, kappa_c = kappa_c, kappa_l = kappa_l)

	# guess value functions
	guess_dfrm = data.frame(grid_dfrm, 
				S_utility=ss_utility(ss_args_s, discount_factor), 
				I_utility=ss_utility(ss_args_i, discount_factor)*phi, 
				R_utility=ss_utility(ss_args_r, discount_factor))
    ## initialize others_choices
    others_choices_init = data.frame(S=rep(l_S_guess,length.out=nrow(grid_dfrm)), I=rep(l_I_guess,length.out=nrow(grid_dfrm)), R=rep(l_R_guess,length.out=nrow(grid_dfrm)))
    others_choices = others_choices_init

    SWF_guess = rowSums(guess_dfrm[,1:3]*guess_dfrm[,4:6])
}

if(hot_start == 1 & length(scriptargs) != 0) {
	message("Pulling up hot start guess...")
	list_of_previous_solves = list.files("../../Results/value_policy_functions/")
	
	chosen_guess = list_of_previous_solves[which(list_of_previous_solves==paste0("value_policy_functions__",scenario_label,".csv"))]
	# message("Nearest file to start from is ",chosen_guess)

	previous_solve = read.csv(paste0("../../Results/value_policy_functions/",chosen_guess))

	if(nrow(previous_solve) != nrow(grid_dfrm)) {
		previous_solve = grid_grow(previous_solve, grid_dfrm, simplex_check)
	}

	u_S_guess = previous_solve$lifetime_utility_S
	u_I_guess = previous_solve$lifetime_utility_I
	u_R_guess = previous_solve$lifetime_utility_R

	l_S_guess = previous_solve$labor_supply_S
	l_I_guess = previous_solve$labor_supply_I
	l_R_guess = previous_solve$labor_supply_R

	guess_dfrm = data.frame(grid_dfrm, 
					S_utility=u_S_guess, 
					I_utility=u_I_guess, 
					R_utility=u_R_guess)

	SWF_guess = rowSums(guess_dfrm[,1:3]*guess_dfrm[,4:6])

    ## initialize others_choices
    others_choices_init = data.frame(S=l_S_guess, I=l_I_guess, R=l_R_guess)
    others_choices = others_choices_init
}

contval_list_init = list(S=guess_dfrm$S_utility, I=guess_dfrm$I_utility, R=guess_dfrm$R_utility, SWF=SWF_guess)
contval_list = contval_list_init

#####
# D: run VFI solver loop. layout inside "dp_solver()" is:
# outermost loop (while): value function iteration. updates the continuation values with results from inner loops
# middle loop (foreach): iterates over states. this loop is parallelized since problems at each state are independent of each other state. middle loop is split into batches to facilitate convergence.
# innermost loop: consistency of choices. this loop makes sure all types' choices are consistent with "best-responding" to each other's choices
#####

policy_list = list()
policy_list$switch = policy #1 if policy will be put in place
policy_list$start_date = 25 #policy start time in days
policy_list$end_date = 31 #policy end time in days
policy_list$labor_supply_level = 0.7 #reduce labor supply to this much of initial (day 1) level
policy_list$start_trigger = 0.015 #start 1-day stay-at-home if infection crosses this threshold

total_solve_time = proc.time()[3]
solved_values = dp_solver()		# the actual solver call!
total_solve_time = proc.time()[3] - total_solve_time
message("Total solver run time: ", round(total_solve_time/60,3), " minutes")
solved_values = cbind(solved_values, grid_dfrm, gamma_c=gamma_c, gamma_l=gamma_l, risk_aversion=risk_aversion, s=s,alpha_U=alpha_U,phi=phi,rho_c=rho_c,rho_l=rho_l,rho_o=rho_o,duration=duration,policy=0)

contval_list = list(S=solved_values$lifetime_utility_S, I=solved_values$lifetime_utility_I, R=solved_values$lifetime_utility_R, SWF=solved_values$SWF)

# save output as a file
fwrite(solved_values, file=paste0("../../Results/value_policy_functions/value_policy_functions__",scenario_label,".csv"))

#####
# C: generate baseline time series
####
SIR_baseline_solve = solved_values
SIR_baseline_solve$labor_supply_S = rep(Lbm,length.out=nrow(solved_values))
SIR_baseline_solve$labor_supply_I = rep(Lbm,length.out=nrow(solved_values))
SIR_baseline_solve$labor_supply_R = rep(Lbm,length.out=nrow(solved_values))
SIR_baseline_solve$consumption_S = rep(wage*Lbm,length.out=nrow(solved_values))
SIR_baseline_solve$consumption_I = rep(wage*phi*Lbm,length.out=nrow(solved_values))
SIR_baseline_solve$consumption_R = rep(wage*Lbm,length.out=nrow(solved_values))

tps_time = proc.time()[3]
results_just_SIR = generate_time_series(SIR_baseline_solve, interpolator="multilinear")
tps_time = proc.time()[3] - tps_time
message("Time to generate pure SIR multilinear interpolation: ", round(tps_time/60,3), " minutes")
solved_values = cbind(solved_values, grid_dfrm, gamma_c=gamma_c, gamma_l=gamma_l, risk_aversion=risk_aversion, s=s,policy=0)

tps_time = proc.time()[3]
results_nopolicy_linear = generate_time_series(solved_values, interpolator="multilinear")
tps_time = proc.time()[3] - tps_time
message("Time to generate behavioral multilinear interpolation: ", round(tps_time/60,3), " minutes")

if(policy_list$switch==1){
	results_SIR_policy = generate_time_series_policy(SIR_baseline_solve,policy_list)
	results_policy = generate_time_series_policy(solved_values,policy_list)
}

# stopCluster(optim_cl)
stopCluster(cl)

message("In the no-policy pure-SIR scenario, infection prevalence peaks in day ", which.max(results_just_SIR$I), " at ", 100*round(max(results_just_SIR$I),3), "% of the population.")

# message("In the no-policy behavioral (tps) scenario, infection prevalence peaks in day ", which.max(results_nopolicy$I), " at ", 100*round(max(results_nopolicy$I),3), "% of the population.")
# fwrite(results_nopolicy, file=paste0("../../Results/time_series/csvs/SIR_series__contact_u_",gamma_c,"_",gamma_l,"_crra_",risk_aversion,"_tps.csv"))

message("In the no-policy behavioral (linear) scenario, infection prevalence peaks in day ", which.max(results_nopolicy_linear$I), " at ", 100*round(max(results_nopolicy_linear$I),3), "% of the population.")
fwrite(results_nopolicy_linear, file=paste0("../../Results/time_series/csvs/",scenario_label,"_",gamma_c,"_",gamma_l,"_crra_",risk_aversion,"_gridsize_",S_gridlength,".csv"))

fwrite(results_just_SIR, file=paste0("../../Results/time_series/csvs/SIR_series__just_SIR.csv"))

#####
# D: plot projected series
#####

ts_plotlist_nopolicy_linear = generate_plots(results_nopolicy_linear)
ts_plotlist_justSIR = generate_plots(results_just_SIR)

png(paste0("../../Results/time_series/plots/justSIR_results.png"), width=1600, height=900)
plot_grid(plotlist=ts_plotlist_justSIR, nrow=3)
dev.off()
png(paste0("../../Results/time_series/plots/",scenario_label,"_",I_gridlength,".png"), width=1600, height=900)
plot_grid(plotlist=ts_plotlist_nopolicy_linear, nrow=3)
dev.off()

if(policy_list$switch==0){
	fwrite(results_nopolicy_linear, file=paste0("../../Results/Figure2_Models_epi/",scenario_label,".csv"))
}

if(policy_list$switch==1){
	message("In the behavioral scenario with policy, infection prevalence peaks in day ", which.max(results_policy$I), " at ", 100*round(max(results_policy$I),3), "% of the population.")

	fwrite(results_policy, file=paste0("../../Results/time_series/csvs/policy_SIR_series__contact_u_",gamma_c,"_",gamma_l,"_crra_",risk_aversion,"_EZ_",ez_prefs,"_rho_",rho,".csv"))

	fwrite(cbind(results_policy,policy=1), file=paste0("../../Results/Figure3_Policy_scenarios/",scenario_label,"_policy.csv"))
	fwrite(cbind(results_nopolicy,policy=0), file=paste0("../../Results/Figure3_Policy_scenarios/",scenario_label,"_nopolicy.csv"))

	fwrite(cbind(results_SIR_policy,policy=1), file=paste0("../../Results/Figure3_Policy_scenarios/pure_SIR_policy.csv"))
	fwrite(cbind(results_just_SIR,policy=0), file=paste0("../../Results/Figure3_Policy_scenarios/pure_SIR_nopolicy.csv"))


	ts_plotlist_policy = generate_plots(results_policy)
	png(paste0("../../Results/time_series/plots/policy_",scenario_label,".png"), width=1600, height=900)
	plot_grid(plotlist=ts_plotlist_policy, nrow=3)
	dev.off()
}

solved_values$simplex_check = ifelse(solved_values$S + solved_values$I + solved_values$R<=1, 1, 0)

png(paste0("../../Results/value_policy_functions/values_and_policies_",scenario_label,"_",I_gridlength,".png"), width=1200, height=800)
par(mfrow=c(2,3))
chart_value_fn(df=solved_values, type="S", var="lifetime_utility")
chart_value_fn(df=solved_values, type="I", var="lifetime_utility")
chart_value_fn(df=solved_values, type="R", var="lifetime_utility")
chart_value_fn(df=solved_values, type="S", var="labor_supply")
chart_value_fn(df=solved_values, type="I", var="labor_supply")
chart_value_fn(df=solved_values, type="R", var="labor_supply")
dev.off()
}