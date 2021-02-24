##############################################################################
# Functions to calibrate epicon model
##############################################################################

## labor supply calculator
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

# labor_supply(parms, benchmarks)
# labor_supply(c(0.3,parms[2],parms[3],parms[4]), benchmarks)
# labor_supply(parms, benchmarks, state=state, type="dynamic")
# labor_supply(c(0.3,parms[2],parms[3],parms[4]), benchmarks, state=state, type="dynamic")

## wage elasticity of labor supply
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

# wels(parms,benchmarks)

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

## price premium closeness
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

## wage premium closeness
wpm_cost <- function(parms, benchmarks) {

	s = as.numeric(parms[1])
	alpha_U = as.numeric(parms[2])
	gamma_c = as.numeric(parms[3])
	gamma_l = as.numeric(parms[4])

	p = benchmarks$p
	Cbm = benchmarks$Cbm
	Lbm = benchmarks$Lbm
	wageprem = benchmarks$wageprem
	w = benchmarks$w
	rho_c = benchmarks$rho_c
	rho_l = benchmarks$rho_l
	risk_aversion = 1

	l_with = labor_supply(parms,benchmarks)

	benchmarks_without = benchmarks
	benchmarks_without$w = benchmarks$w*benchmarks$wageprem
	parms_without = c(s, alpha_U, gamma_c, 0)
	l_without = labor_supply(parms_without,benchmarks_without)

	c_with = w*l_with
	c_without = w*l_without
	
	u_with_contact_and_premium = total_utility(c_with, l_with, Cbm, Lbm, gamma_c, gamma_l, rho_c, rho_l, risk_aversion, s, alpha_U)
	u_without_contact_or_premium = total_utility(c_without, l_without, Cbm, Lbm, gamma_c, 0, rho_c, rho_l, risk_aversion, s, alpha_U)

	cost = (u_with_contact_and_premium - u_without_contact_or_premium)^2

	return(cost)
}

calibration_objective <- function(weights,parms,benchmarks) {

	weighted_cost = weights[1]*ls_cost(parms, benchmarks) +
					weights[2]*elw_cost(parms, benchmarks) +
					weights[3]*ppm_cost(parms, benchmarks) +
					weights[4]*wpm_cost(parms, benchmarks) +
					weights[5]*budget_cost(parms, benchmarks)

	return(weighted_cost)
}

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
