##############################################################################
# Script to calibrate structural utility parameters
##############################################################################

calibration_time = proc.time()[3]

# Calibration search hyperparameters
# Some known "good" grid points:
# Lbar = 24, income = 58000: s_grid = seq(1.221, 1.229, length.out=200), alpha_U_grid = seq(0.199, 0.209, length.out=200). BC gap = 0.08640008.
# Lbar = 12, income = 58000: s_grid = seq(1.1, 1.6, length.out=200), alpha_U_grid = seq(0.1, 0.5, length.out=200). BC gap = 0.0000007595356
# Lbar = 12, income = 44000: s_grid = seq(1.1, 1.6, length.out=200), alpha_U_grid = seq(0.1, 0.5, length.out=200). BC gap = 0.000001173954

s_grid = seq(1.1, 1.8, length.out=200)
alpha_U_grid = seq(0.1, 0.5, length.out=200)
gamma_c_grid = 0
gamma_l_grid = 0
wageprem = 1
priceprem = 1
elasticity = 0.15

weights = c(1,1,0,0,1) # weights on calibration targets: labor supply, elasticity of labor supply, price premium, wage premium, budget constraint

if(social==1) {
	weights = c(1,1,1,1,1)

	wageprem = 1.1
	priceprem = 0.9

	s_grid = seq(1.1, 1.6, length.out=50)
	alpha_U_grid = seq(0.1, 0.5, length.out=50)
	gamma_c_grid = seq(from=0, to=0.006, length.out=75)
	gamma_l_grid = seq(from=0, to=0.011, length.out=75)
}

uniroot_lower_scaler = 100
uniroot_upper_scaler = 10

# Givens/benchmarks
p = 1
benchmarks = data.frame(risk_aversion = risk_aversion, tau_parm=tau_parm, discount_factor=discount_factor, rho_c=rho_c, rho_l=rho_l, rho_o=rho_o, p=p, Cbm=Cbm, Lbm=Lbm, w=wage, wageprem = wageprem, priceprem = priceprem, elasticity = elasticity, pi_r = pi_r, pi_d = pi_d, Lbar = Lbar, nonlabor_income = nonlabor_income, phi = phi)

write.csv(benchmarks, file="calibration_benchmarks.csv")

# VSL inputs/targets
delta_pi_d = 1/10000		# increase in probability of death
delta_c = VSL*delta_pi_d	# compensation for increased death risk

#####
# Period utility calibration
#####

##### Do grid search for parameter values (search_again==1) or load pre-computed values (search_again==0)

if(calibrate==0) {
	best_parms = read.csv("calibrated_parameters.csv")[,-1]
	values_at_best_parms = 
		data.frame(labor_supply = labor_supply(best_parms,benchmarks)) %>%
		mutate(labor_supply_I = labor_supply(best_parms,benchmarks,type="I")) %>%
		mutate(consumption = labor_supply*best_parms$w) %>%
		mutate(labor_supply_elasticity = wels(best_parms,benchmarks)) %>%
		mutate(price_premium_gap = ppm_cost(as.numeric(best_parms),benchmarks)) %>%
		mutate(wage_premium_gap = wpm_cost(as.numeric(best_parms),benchmarks)) %>%
		mutate(budget_constraint_gap = budget_cost(best_parms,benchmarks)) %>%
		mutate(risk_aversion = risk_aversion)
}

if(calibrate==1) {
	parms_grid = expand.grid(s=s_grid, alpha_U=alpha_U_grid, gamma_c=gamma_c_grid, gamma_l=gamma_l_grid)
	nrow(parms_grid)
	print(benchmarks)
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

		for(k in 1:nrow(split_grid[[i]]) ) {
			
			parms_vector = split_grid[[i]][k,]
			split_result[k] = calibration_objective(weights=weights, parms=parms_vector, benchmarks=benchmarks)

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
		labor_supply_I = labor_supply(best_parms,benchmarks,type="I"),
		consumption=labor_supply(best_parms,benchmarks)*best_parms$w,
		labor_supply_elasticity=wels(best_parms,benchmarks),
		price_premium_gap=ppm_cost(as.numeric(best_parms),benchmarks),
		wage_premium_gap=wpm_cost(as.numeric(best_parms),benchmarks),
		budget_constraint_gap=budget_cost(best_parms,benchmarks),
		loss=min(loss),
		risk_aversion = risk_aversion)

	if(best_parms$gamma_c==0&best_parms$gamma_l==0) {
		write.csv(best_parms, file="calibrated_parameters_nosocial.csv")
		write.csv(values_at_best_parms, file="targets_at_calibrated_parameters_nosocial.csv")
	}
}

message("Best parameters: ")
print(best_parms)
message("Values at best parameters: ")
print(values_at_best_parms)

#####
# VSL calibration
#####

benchmarks_I = benchmarks
benchmarks_I$w = wage*phi
l_I = labor_supply(best_parms, benchmarks_I)
c_I = wage*phi*l_I

s = best_parms$s
alpha_U = best_parms$alpha_U

U_R = total_utility(Cbm,Lbm,Cbm,Lbm,gamma_c=0,gamma_l=0,rho_c=0,rho_l=0,risk_aversion=risk_aversion,s=s,alpha=alpha_U, Lbar=Lbar)/(1-discount_factor)

total_utility(Cbm,Lbm,Cbm,Lbm,gamma_c=0,gamma_l=0,rho_c=0,rho_l=0,risk_aversion=risk_aversion,s=s,alpha=alpha_U, Lbar=24)/(1-discount_factor)
total_utility(Cbm,Lbm,Cbm,Lbm,gamma_c=0,gamma_l=0,rho_c=0,rho_l=0,risk_aversion=risk_aversion,s=s,alpha=alpha_U, Lbar=12)/(1-discount_factor)

message("Lower VSL_balance endpoint value")
print(VSL_balance(-uniroot_lower_scaler*U_R, c_I, pi_d, pi_r, discount_factor, U_R, risk_aversion, s, alpha_U, C=0, L=0, delta_c = delta_c, delta_pi_d = delta_pi_d, Lbar = Lbar))

message("Upper VSL_balance endpoint value")
print(VSL_balance(uniroot_upper_scaler*U_R, c_I, pi_d, pi_r, discount_factor, U_R, risk_aversion, s, alpha_U, C=0, L=0, delta_c = delta_c, delta_pi_d = delta_pi_d, Lbar = Lbar))

(result <- uniroot(VSL_balance, interval=c(-uniroot_lower_scaler*U_R,uniroot_upper_scaler*U_R), c_I = c_I, pi_d = pi_d, pi_r = pi_r, discount_factor = discount_factor, U_R = U_R, risk_aversion = risk_aversion, s = s, alpha_U = alpha_U, C=0, L=0, delta_c = delta_c, delta_pi_d = delta_pi_d, Lbar = Lbar))

U_D = result$root

# implied_VSL_balance(1e10, c_I = c_I, pi_d = pi_d, pi_r = pi_r, discount_factor = discount_factor, U_R = U_R, risk_aversion = 1.1, s = s, alpha_U = alpha_U, C=0, L=0, U_D = 0, delta_pi_d = delta_pi_d, Lbar = Lbar)

# U_R = total_utility(Cbm,Lbm,Cbm,Lbm,gamma_c=0,gamma_l=0,rho_c=0,rho_l=0,risk_aversion=1.2,s=s,alpha=alpha_U, Lbar=Lbar)/(1-discount_factor)

# (implied_VSL_result <- uniroot(implied_VSL_balance, interval=c(0, 1e50), c_I = c_I, pi_d = pi_d, pi_r = pi_r, discount_factor = discount_factor, U_R = U_R, risk_aversion = 1.2, s = s, alpha_U = alpha_U, C=0, L=0, U_D = 0, delta_pi_d = delta_pi_d, Lbar = Lbar))

# (implied_VSL = implied_VSL_result$root/delta_pi_d)

message("Time endowment is ", Lbar, " hours.")

message("Lifetime utility of death (U_D): ")
print(U_D)

U_I_value = U_I(result$root, c_I, pi_d, pi_r, discount_factor, U_R, risk_aversion, s, alpha_U, C=Cbm, L=Lbm, Lbar)
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

calibration_time = round((proc.time()[3] - calibration_time)/60,3)

message("Total calibration time: ", calibration_time, " minutes.")
