##### Script to generate U_D to calibrate model to hit target VSL
# This script targets the I-type's VSL by setting U_D. Targeting the I-type's VSL has some advantages over targeting the S-type's VSL: since the I-type directly faces the risk of death, this allows us to calibrate U_D without requiring any P_I calculations.
# This is only designed to work for the base (no social) model, but for any level of risk aversion

rm(list=ls())

library(ggplot2)
library(cowplot)
library(doParallel)
library(chebpol)
library(pracma)
library(rootSolve)
library(compiler)
library(data.table)
library(dplyr)
library(fields)
library(progress)

options("scipen"=100)

enableJIT(3) 

source("functions.R")

Cbm = (58000/365)	# daily annuitized consumption flow in dollars
Lbm = 0.3333*24		# hours of day spent working
A = Cbm/Lbm 								# hourly wage

VSL = 10000000
delta_pi_d = 1/10000
delta_c = VSL*delta_pi_d

discount_rate = 0.04								# discount rate
discount_factor = (1/(1+discount_rate))^(1/365)				# dayly discount factor

duration = 5.1						# duration of infectious period (rename duration as duration)
est_cfr = 0.015					# conditional on leaving I, 5% chance end up in D; Pr(I->D)
proportion_removed = 1/duration		#1 - (1 - (1/duration))^7
pi_r = (1-est_cfr)*proportion_removed
pi_d = est_cfr*proportion_removed

alpha_U = 0.3684211					# CES share parameter
gamma_c = 0								# consumption contact utility parameter
gamma_l = 0								# labor contact utility parameter
s =  1.457895 								# elasticity of substitution between leisure and consumption
Lbar_endow = 12

eta = 0.1

phi = 0.8554632

l_I = 7.80018421836839							# taken from solved policy functions -- ideally should be taken from solving the I-type's static optimization problem
c_I = A*phi*l_I

###### Utilities

U_R = total_utility(Cbm,Lbm,Cbm,Lbm,gamma_c=0,gamma_l=0,rho_c=0,rho_l=0,eta=eta,s=s,alpha=alpha_U, Lbar=Lbar_endow)/(1-discount_factor)

U_I <- function(U_D, c_I, pi_d, pi_r, discount_factor, U_R, eta, s, alpha_U, Cbm=0, Lbm=0, delta_c=0, delta_pi_d=0, Lbar=Lbar_endow) {
	u_it = total_utility(c_I+delta_c, l_I,Cbm,Lbm,gamma_c=0,gamma_l=0,rho_c=0,rho_l=0,eta=eta,s=s,alpha=alpha_U, Lbar=Lbar_endow)

	U_I_normal = (total_utility(c_I, l_I,Cbm,Lbm,gamma_c=0,gamma_l=0,rho_c=0,rho_l=0,eta=eta,s=s,alpha=alpha_U, Lbar=Lbar) + discount_factor*(pi_r*U_R + pi_d*U_D))/(1 - discount_factor*(1 - pi_r - pi_d))
	
	value = u_it + discount_factor*(pi_r*U_R + (pi_d+delta_pi_d)*U_D + (1 - pi_r - (pi_d+delta_pi_d))*U_I_normal)
	
	return(value)
}

VSL_balance <- function(U_D, c_I, pi_d, pi_r, discount_factor, U_R, eta, s, alpha_U, Cbm=0, Lbm=0, delta_c, delta_pi_d, Lbar=Lbar_endow) {
	balance = U_I(U_D, c_I, pi_d, pi_r, discount_factor, U_R, eta, s, alpha_U, Cbm=0, Lbm=0, delta_c = 0, delta_pi_d = 0) - U_I(U_D, c_I, pi_d, pi_r, discount_factor, U_R, eta, s, alpha_U, Cbm=0, Lbm=0, delta_c = delta_c, delta_pi_d = delta_pi_d, Lbar = Lbar_endow)
	return(balance)
}

VSL_balance(0, c_I, pi_d, pi_r, discount_factor, U_R, eta, s, alpha_U, Cbm=0, Lbm=0, delta_c = delta_c, delta_pi_d = delta_pi_d, Lbar = Lbar_endow)

VSL_balance(-2883300.45, c_I, pi_d, pi_r, discount_factor, U_R, eta, s, alpha_U, Cbm=0, Lbm=0, delta_c = delta_c, delta_pi_d = delta_pi_d, Lbar = Lbar_endow)

(result <- uniroot(VSL_balance, interval=c(-U_R*10,U_R), c_I = c_I, pi_d = pi_d, pi_r = pi_r, discount_factor = discount_factor, U_R = U_R, eta = eta, s = s, alpha_U = alpha_U, Cbm=0, Lbm=0, delta_c = delta_c, delta_pi_d = delta_pi_d, Lbar = Lbar_endow))

print(result$root)
print(U_I(result$root, c_I, pi_d, pi_r, discount_factor, U_R, eta, s, alpha_U, Cbm=0, Lbm=0))
print(U_R)
