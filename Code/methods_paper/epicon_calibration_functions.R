##############################################################################
# Functions to calibrate epicon model
##############################################################################

# Analytical formulas for calibration -- taken from Dan's notes on CEs function calibration. NOTE that alpha_hat = leisure intensity parameter, alpha = leisure share parameters. alpha_hat = (1-alpha)/alpha.

## analytical formula for labor supply
l_star <- function(sigma,alpha_hat,wage,Lbar) {
	output <- Lbar - alpha_hat*Lbar*wage/(wage^sigma + alpha_hat*wage)
	return(output)
}

## analytical formula for uncompensated elasticity of labor supply wrt wage
eta_lw <- function(sigma,alpha_hat,wage,Lbar) {
	output <- alpha_hat*wage*(sigma-1)/(wage^sigma + wage*alpha_hat) 
	return(output)
}

## analytical formula calibration objective
calibration_objective_analytical <- function(calibration_parms,calibration_targets,extra_parms) {

	sigma <- calibration_parms[1]
	alpha <- calibration_parms[2]
	alpha_hat <- (1-alpha)/alpha

	target_laborsupply <- calibration_targets$labor_supply
	target_elasticity <- calibration_targets$elasticity

	wage <- extra_parms$wage
	Lbar <- extra_parms$Lbar

	weighted_cost = 1e6*(l_star(sigma,alpha_hat,wage,Lbar) - target_laborsupply)^2 +
					1e6*(eta_lw(sigma,alpha_hat,wage,Lbar) - target_elasticity)^2

	return(weighted_cost)
}

## labor supply calculator. Calculates labor supply which maximizes whatever utility is given to optim. "TU_wrapper" is spelled out in functions.R, and is the static per-period utility function with no transitions. This is meant to reflect the pre-epidemic steady state, when there's only 1 type and no transitions to any other type-state. (It's relatively straightforward to show that the lifetime utility function under CES preferences is itself a CES function, so you could justify this procedure as maximizing lifetime utility while ignoring a multiplicative constant.)
labor_supply <- function(parms, benchmarks, type="S", model="decentralized") {

	rho_c = benchmarks$rho_c
	rho_l = benchmarks$rho_l

	C = benchmarks$Cbm
	L = benchmarks$Lbm
	w = benchmarks$w
	Lbar = benchmarks$Lbar
	nonlabor_income = benchmarks$nonlabor_income
	phi = benchmarks$phi

	s = parms[1]
	alpha_U = parms[2]
	gamma_c = parms[3]
	gamma_l = parms[4]

	if(type=="S") {phi = 1}

	result = optim(par = Lbm, fn = TU_wrapper, C = C, L = L, gamma_c = gamma_c, gamma_l = gamma_l, rho_c = rho_c, rho_l = rho_l, risk_aversion = 1, s = s, alpha = alpha_U, w=w, Lbar=Lbar, nonlabor_income = nonlabor_income, phi = phi, control = list(fnscale=-1), method = "L-BFGS-B", lower=0.01, upper=Lbar)
	labor = result$par
	return(labor)
}

## wage elasticity of labor supply -- calculated numerically rather than by analytical formula. Might be a good idea to derive the formula and plug it in here for clarity over whether we're picking up compensated or uncompensated.
wels <- function(parms, benchmarks) {
	l = labor_supply(parms, benchmarks)
	w = benchmarks$w

	benchmarks_p = benchmarks
	perturbation = 0.1
	benchmarks_p$w = benchmarks_p$w + perturbation
	l_perturbed = labor_supply(parms,benchmarks_p)

	dldw = (l_perturbed - l)/perturbation

	elw = (w/l)*dldw
	return(elw)
}

## labor supply closeness
ls_cost <- function(parms, benchmarks) {
	calc_ls = labor_supply(parms,benchmarks)
	target_ls = benchmarks$Lbm

	cost = (calc_ls - target_ls)^2

	return(cost)
}

## budget constraint closeness
budget_cost <- function(parms,benchmarks) {
	calc_ls = labor_supply(parms,benchmarks)
	budget_gap = calc_ls*benchmarks$w - benchmarks$Cbm

	cost = budget_gap^2

	return(cost)
}

## elasticity closeness
elw_cost <- function(parms, benchmarks) {
	calc_elw = wels(parms, benchmarks)*100
	target_elw = benchmarks$elasticity*100

	cost = (calc_elw - target_elw)^2

	return(cost)
}

## price premium closeness. Only needed for social utility.
ppm_cost <- function(parms, benchmarks) {

	s = as.numeric(parms[1])
	alpha_U = as.numeric(parms[2])
	gamma_c = as.numeric(parms[3])
	gamma_l = as.numeric(parms[4])

	p = benchmarks$p
	Cbm = benchmarks$Cbm
	Lbm = benchmarks$Lbm
	priceprem = benchmarks$priceprem
	w = benchmarks$w
	rho_c = benchmarks$rho_c
	rho_l = benchmarks$rho_l
	risk_aversion = 1

	l_with = labor_supply(parms,benchmarks)

	benchmarks_without = benchmarks
	benchmarks_without$p = benchmarks$p*benchmarks$priceprem
	parms_without = c(s, alpha_U, 0, gamma_l)
	l_without = labor_supply(parms_without,benchmarks_without)

	c_with = w*l_with
	c_without = w*l_without

	u_with_contact_and_premium = total_utility(c_with, l_with, Cbm, Lbm, gamma_c, gamma_l, rho_c, rho_l, risk_aversion, s, alpha_U)
	u_without_contact_or_premium = total_utility(c_without, l_without, Cbm, Lbm, 0, gamma_l, rho_c, rho_l, risk_aversion, s, alpha_U)

	cost = (u_with_contact_and_premium - u_without_contact_or_premium)^2

	return(cost)
}

## wage premium closeness. Only needed for social utility.
wpm_cost <- function(parms, benchmarks) {

	s = as.numeric(parms[1])
	alpha_U = as.numeric(parms[2])
	gamma_c = as.numeric(parms[3])
	gamma_l = as.numeric(parms[4])

	p = benchmarks$p
	Cbm = benchmarks$Cbm
	Lbm = benchmarks$Lbm
	wageprem = 1#benchmarks$wageprem
	w = benchmarks$w
	rho_c = benchmarks$rho_c
	rho_l = benchmarks$rho_l
	risk_aversion = 1

	l_with = labor_supply(parms,benchmarks)

	benchmarks_without = benchmarks
	# benchmarks_without$w = benchmarks$w*benchmarks$wageprem # disabled on 22/05/2021 because stupid script kept saying "wageprem not found" EVEN WHEN IT WAS TOTALLY SET TO 1
	benchmarks_without$w = benchmarks$w
	parms_without = c(s, alpha_U, gamma_c, 0)
	l_without = labor_supply(parms_without,benchmarks_without)

	c_with = w*l_with
	c_without = w*l_without
	
	u_with_contact_and_premium = total_utility(c_with, l_with, Cbm, Lbm, gamma_c, gamma_l, rho_c, rho_l, risk_aversion, s, alpha_U)
	u_without_contact_or_premium = total_utility(c_without, l_without, Cbm, Lbm, gamma_c, 0, rho_c, rho_l, risk_aversion, s, alpha_U)

	cost = (u_with_contact_and_premium - u_without_contact_or_premium)^2

	return(cost)
}

# The actual objective function for calibration. We're treating the calibration process as a weighted "cost" minimization problem, using "cost" in the engineering/optimization sense (i.e. as a "loss" or "error" function) rather than the economic one (an actual cost). Weights are user-supplied. 
calibration_objective <- function(weights,parms,benchmarks) {

	weighted_cost = weights[1]*ls_cost(parms, benchmarks) +
					weights[2]*elw_cost(parms, benchmarks) +
					weights[3]*ppm_cost(parms, benchmarks) +
					weights[4]*wpm_cost(parms, benchmarks) +
					weights[5]*budget_cost(parms, benchmarks)

	return(weighted_cost)
}

# This was something we sketched out way back when as an option to calibrate risk aversion. We decided against it since we didn't want to bake this in.
infection_elasticity_of_labor_supply <- function(parms,benchmarks,state) {
	state_new = state
	state_new$I = state_new$I*1.01
	
	dynls = labor_supply(parms, benchmarks, state, type="dynamic")
	dynls_new = labor_supply(parms, benchmarks, state_new, type="dynamic")

	dl_dI = (dynls - dynls_new)/(state_new$I - state$I)

	elasticity = (state_new$I/dynls)*dl_dI

	return(elasticity)
}
# infection_elasticity_of_labor_supply(parms,benchmarks,state)
# infection_elasticity_of_labor_supply(c(0.3,parms[2],parms[3],parms[4]),benchmarks,state)

U_I <- function(U_D, c_I, pi_d, pi_r, discount_factor, U_R, risk_aversion, s, alpha_U, C=0, L=0, delta_c=0, delta_pi_d=0, Lbar=24) {
	u_it = total_utility(c_I+delta_c, l_I, C, L, gamma_c=0,gamma_l=0,rho_c=0,rho_l=0,risk_aversion=risk_aversion,s=s,alpha=alpha_U,Lbar=Lbar) # the one-period value of getting that extra consumption

	modified_discount_factor = (1 - discount_factor*(1 - pi_r - pi_d)) # delta_pi_d changes for just one period, so this uses the standard pi_d instead of new_pi_d
	new_pi_d = pi_d + delta_pi_d # this is the one-period-higher chance of death

	U_I_normal = (total_utility(c_I, l_I, C, L, gamma_c=0,gamma_l=0,rho_c=0,rho_l=0,risk_aversion=risk_aversion,s=s,alpha=alpha_U, Lbar=Lbar) + discount_factor*(pi_r*U_R + pi_d*U_D))/modified_discount_factor # this is the lifetime utility of being infected with the usual probabilities and consumption
	
	value = u_it + discount_factor*(pi_r*U_R + new_pi_d*U_D + (1 - pi_r - new_pi_d)*U_I_normal) # the lifetime value of one-period higher consumption and higher deathrisk
	
	return(value)
}

VSL_balance <- function(U_D, c_I, pi_d, pi_r, discount_factor, U_R, risk_aversion, s, alpha_U, C=0, L=0, delta_c, delta_pi_d, Lbar) {
	balance = U_I(U_D, c_I, pi_d, pi_r, discount_factor, U_R, risk_aversion, s, alpha_U, C=C, L=L, delta_c = 0, delta_pi_d = 0, Lbar = Lbar) - U_I(U_D, c_I, pi_d, pi_r, discount_factor, U_R, risk_aversion, s, alpha_U, C=C, L=L, delta_c = delta_c, delta_pi_d = delta_pi_d, Lbar = Lbar)
	return(balance)
}

implied_VSL_balance <- function(delta_c, U_D, c_I, pi_d, pi_r, discount_factor, U_R, risk_aversion, s, alpha_U, C=0, L=0,  delta_pi_d, Lbar) {
	balance = U_I(U_D, c_I, pi_d, pi_r, discount_factor, U_R, risk_aversion, s, alpha_U, C=C, L=L, delta_c = 0, delta_pi_d = 0, Lbar = Lbar) - U_I(U_D, c_I, pi_d, pi_r, discount_factor, U_R, risk_aversion, s, alpha_U, C=C, L=L, delta_c = delta_c, delta_pi_d = delta_pi_d, Lbar = Lbar)
	return(balance)
}
