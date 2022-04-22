##############################################################################
# Script to calibrate epicon model. Modeled off the Mathematica notebook "contact_utility_calibration_2.nb"
##############################################################################

rm(list=ls())

library(doSNOW)
library(ggplot2)
library(cowplot)
source("functions.R")
source("epicon_calibration_functions.R")

scriptargs = commandArgs(trailingOnly = TRUE)

search_again = 1
ncores = 80
social = 0

# Calibration search hyperparameters
# Some known "good" grid points:
# Lbar = 24, income = 58000: s_grid = seq(1.221, 1.229, length.out=200), alpha_U_grid = seq(0.199, 0.209, length.out=200). BC gap = 0.08640008.
# Lbar = 12, income = 58000: s_grid = seq(1.1, 1.6, length.out=200), alpha_U_grid = seq(0.1, 0.5, length.out=200). BC gap = 0.0000007595356
# Lbar = 12, income = 44000: s_grid = seq(1.1, 1.6, length.out=200), alpha_U_grid = seq(0.1, 0.5, length.out=200). BC gap = 0.000001173954
s_grid = seq(1.1, 1.6, length.out=300)
alpha_U_grid = seq(0.1, 0.5, length.out=300)
weights = c(1,1,0,0,1) # weights on calibration targets: labor supply, elasticity of labor supply, price premium, wage premium, budget constraint
wageprem = 1
priceprem = 1
elasticity = 0.15

gamma_c_grid = 0
gamma_l_grid = 0
if(social==1) {
	weights = c(1,1,1,1,1)

	# Targets (excluding Lbm since it's defined above)
	wageprem = 1.1
	priceprem = 0.9
	s_grid = seq(1.1, 1.6, length.out=40)
	alpha_U_grid = seq(0.1, 0.5, length.out=40)
	gamma_c_grid = seq(from=0, to=0.006, length.out=75)
	gamma_l_grid = seq(from=0, to=0.011, length.out=75)
}

uniroot_lower_scaler = 100
uniroot_upper_scaler = 10

# Key utility sensitivity parameters
annual_income = 58000				# annual income (consumption target), in dollars
time_endowment = 12					# endowment of time usable for consumption and leisure, in hours
daily_cons_contacts = 5.166426 		# number of contacts at a consumption activity per day
daily_work_contacts = 7.513051		# number of contacts at a labor activity per day
daily_other_contacts = 3.548544 #as.numeric(scriptargs[11])
eta = 0.1

# Givens/benchmarks
p = 1
discount_rate = 0.04
discount_factor = (1/(1+discount_rate))^(1/365)
Cbm = (annual_income/365)	# daily annuitized consumption flow in dollars
Lbm = 0.3333*24		# hours of day spent working
w = Cbm/Lbm			# the wage comes directly from the budget constraint at the benchmark
phi = 0.8554632		# productivity of infected
nonlabor_income = 0	# money in the bank

rho_c = daily_cons_contacts/(Cbm^2)	# daily_cons_contacts = rho_c*Cbm*Cbm
rho_l = daily_work_contacts/(Lbm^2)	# daily_work_contacts = rho_l*Lbm*Lbm
rho_o = daily_other_contacts		# O*O := 1 (normalization)
duration = 5.1						# duration of infectious period (rename duration as duration)
est_cfr = 0.015					# conditional on leaving I, 5% chance end up in D; Pr(I->D)
proportion_removed = 1/duration		#1 - (1 - (1/duration))^7
pi_r = (1-est_cfr)*proportion_removed
pi_d = est_cfr*proportion_removed
target_R0 = 2.6
tau_parm = target_R0/(duration*(daily_cons_contacts+daily_work_contacts+daily_other_contacts))			# This should be equivalent to just the tau from the compute.beta function in tau.R

# VSL inputs/targets
VSL = 10000000
delta_pi_d = 1/10000
delta_c = VSL*delta_pi_d

# SIR initial state
I_0 = 0
state = data.frame(S=1-I_0,I=I_0,R=0)

#####
# Period utility calibration
#####

##### Set up vector of benchmark values and targets

benchmarks = data.frame(eta=eta, tau_parm=tau_parm, discount_factor=discount_factor, rho_c=rho_c, rho_l=rho_l, rho_o=rho_o, p=p, Cbm=Cbm, Lbm=Lbm, w=w, wageprem = wageprem, priceprem = priceprem, elasticity = elasticity, state = state, pi_r = pi_r, pi_d = pi_d, Lbar = time_endowment, nonlabor_income = nonlabor_income, risk_aversion = eta)

write.csv(benchmarks, file="calibration_benchmarks_nosocial.csv")

##### Do grid search for parameter values (search_again==1) or load pre-computed values (search_again==0)

if(search_again==0) {
	best_parms = read.csv("calibrated_parameters.csv")[,-1]
}

if(search_again==1) {
	parms_grid = expand.grid(s=s_grid, alpha_U=alpha_U_grid, gamma_c=gamma_c_grid, gamma_l=gamma_l_grid)
	nrow(parms_grid)
	message(nrow(parms_grid)," grid points total")

	pb <- txtProgressBar(max = nrow(parms_grid), style = 3)
	progress <- function(n) setTxtProgressBar(pb, n)
	opts <- list(progress = progress)


	bins  <- sort(rep(1:ncores, length.out=nrow(parms_grid)))
	split_grid = split(parms_grid, bins)

	cl = makeCluster(ncores)
	registerDoSNOW(cl)

	loss = foreach(i=seq(1:length(split_grid)), .inorder=TRUE, .options.snow=opts) %dopar% {

		split_result = rep(NA,length.out=nrow(split_grid[[i]]))

		for(k in 1:nrow(split_grid[[i]])) {
			
			parms_vector = as.numeric(split_grid[[i]][k,])
			split_result[k] = calibration_objective(weights, parms_vector, benchmarks)

		}

		split_result
	}

	close(pb)
	stopCluster(cl)

	loss = unlist(loss)

	rownames(loss) = NULL

	parms_and_loss = data.frame(parms_grid, loss=loss)

	best_parms = data.frame(parms_and_loss[which.min(loss),],w=benchmarks$w)
	values_at_best_parms = data.frame(
		labor_supply=labor_supply(best_parms,benchmarks),
		consumption=labor_supply(best_parms,benchmarks)*best_parms$w,
		labor_supply_elasticity=wels(best_parms,benchmarks),
		price_premium_gap=ppm_cost(as.numeric(best_parms),benchmarks),
		wage_premium_gap=wpm_cost(as.numeric(best_parms),benchmarks),
		budget_constraint_gap=budget_cost(best_parms,benchmarks),
		risk_aversion = eta)

	message("Best parameters: ")
	print(best_parms)
	message("Values at best parameters: ")
	print(values_at_best_parms)

	if(best_parms$gamma_c==0&best_parms$gamma_l==0) {
		write.csv(best_parms, file="calibrated_parameters_nosocial.csv")
		write.csv(values_at_best_parms, file="targets_at_calibrated_parameters_nosocial.csv")
	}
}


#####
# VSL calibration
#####

s = best_parms$s
alpha_U = best_parms$alpha_U

U_R = total_utility(Cbm,Lbm,Cbm,Lbm,gamma_c=0,gamma_l=0,rho_c=0,rho_l=0,eta=eta,s=s,alpha=alpha_U, Lbar=time_endowment)/(1-discount_factor)

total_utility(Cbm,Lbm,Cbm,Lbm,gamma_c=0,gamma_l=0,rho_c=0,rho_l=0,eta=eta,s=s,alpha=alpha_U, Lbar=24)/(1-discount_factor)
total_utility(Cbm,Lbm,Cbm,Lbm,gamma_c=0,gamma_l=0,rho_c=0,rho_l=0,eta=eta,s=s,alpha=alpha_U, Lbar=12)/(1-discount_factor)

message("Lower VSL_balance endpoint value")
print(VSL_balance(-uniroot_lower_scaler*U_R, c_I, pi_d, pi_r, discount_factor, U_R, eta, s, alpha_U, C=0, L=0, delta_c = delta_c, delta_pi_d = delta_pi_d, Lbar = time_endowment))

message("Upper VSL_balance endpoint value")
print(VSL_balance(uniroot_upper_scaler*U_R, c_I, pi_d, pi_r, discount_factor, U_R, eta, s, alpha_U, C=0, L=0, delta_c = delta_c, delta_pi_d = delta_pi_d, Lbar = time_endowment))

(result <- uniroot(VSL_balance, interval=c(-uniroot_lower_scaler*U_R,uniroot_upper_scaler*U_R), c_I = c_I, pi_d = pi_d, pi_r = pi_r, discount_factor = discount_factor, U_R = U_R, eta = eta, s = s, alpha_U = alpha_U, C=0, L=0, delta_c = delta_c, delta_pi_d = delta_pi_d, Lbar = time_endowment))

U_D = result$root


message("Time endowment is ", time_endowment, " hours.")

message("Lifetime utility of death (U_D): ")
print(U_D)

U_I_value = U_I(result$root, c_I, pi_d, pi_r, discount_factor, U_R, eta, s, alpha_U, C=Cbm, L=Lbm, time_endowment)
message("Lifetime utility of being infected (U_I): ") 
print(U_I_value)

message("Lifetime utility of being recovered (U_R): ")
print(U_R)

message("Lifetime utility cost of infection: ", round((1 - U_I_value/U_R)*100,3), "%")
message("Lifetime utility cost of death: ", (1 - round(U_D/U_R,3))*100, "%")

best_parms <- data.frame(best_parms, U_I = U_I_value, U_R = U_R)

#####
# Write out calibrated parameters
#####

best_parms$U_D = result$root

write.csv(best_parms, file="calibrated_parameters.csv")
write.csv(values_at_best_parms, file="targets_at_calibrated_parameters.csv")

