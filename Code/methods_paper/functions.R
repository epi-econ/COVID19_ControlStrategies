##############################################################################
# Collection of functions used in epicon script
# Functions contained here:

# simplex_count: given a unit hypercube with n points in each dimension, calculates the number of points which lie in the unit simplex
# Numextract: function to extract numbers from strings
# compute_R0: R0 for a given model run from the initial values (I and choices)
# newly_infected: number of newly infected in period t based on the state and choices
# SIR_update: next-period update for the SIRD system
# tau_t: probability of infection for agent type j based on the state and choices
# private_utility: private utility of agent type j
# contact_utility: contact utility of agent type j
# continuation_value: wrapper for interpolator
# pre_value_fn: unmaximized value function of agent type j
# project_choices: interpolates the choices of each type at the supplied state
# SIR_series: SIR series using projected choices
# make_cheby: takes a uniform grid and converts the nodes to expanded chebyshev
# build_grid: builds a rectangular grid, will be expanded chebyshev nodes if cheby==1, will be a unit hypersimplex instead of a unit hypercube if simplex==1. cheby defaults to 1, simplex defaults to 0 (the latter since ipol() assumes the grid is a hypercube)
# generate_time_series: wrapper for SIR_series
# dp_solver: runs all dynamic programming loops. takes the environment as input ((...) argument)
# generate_plots: wrapper to generate a list of plots from a time series of outcomes
# grid_grow: expands a solved value function using linear interpolation on the smaller grid to project to the larger grid

##############################################################################

simplex_count <- function(n) {
    grid_pieces = build_grid(n, n, n, cheby=0, simplex=0)

    grid_list = grid_pieces[1:3]
    grid_dfrm = grid_pieces$dfrm

    simplex_check = rep(NA,length.out=nrow(grid_dfrm))
    simplex_check = ifelse(rowSums(grid_dfrm)<=1, 1, 0)

    return(sum(simplex_check))
}

Numextract <- function(string){
  as.numeric(unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string))))
}

compute_R0 <- function(parms,choices,SIR) {
    rho_c = parms$rho_c
    rho_l = parms$rho_l
    rho_o = parms$rho_o
    pi_r = parms$pi_r
    pi_d = parms$pi_d
    tau_parm = parms$tau_parm
    phi = parms$phi
    w = parms$w
    kappa_c = parms$kappa_c
    kappa_l = parms$kappa_l

    l_S = choices$l_S
    l_I = choices$l_I
    l_R = choices$l_R
    c_S = choices$l_S*w
    c_I = choices$l_I*w*phi

    S = SIR$S
    I = SIR$I
    R = SIR$R

    parms_to_N_I = data.frame(rho_c = rho_c, rho_l = rho_l, rho_o = rho_o, tau_parm = tau_parm, kappa_c = kappa_c, kappa_l = kappa_l)
    choices_to_N_I = data.frame(c_S = c_S, c_I = c_I, l_S = l_S, l_I = l_I)
    SI_to_N_I = data.frame(S = S, I = I)

    N_I_0 = newly_infected(parms_to_N_I, choices_to_N_I, SI_to_N_I)
    R0 = (N_I_0/I)/(pi_r + pi_d)

    return(R0)
}

# function to calculate the probability type j interacts with an infected
tau_t <- function(type, parms, aggregates) {
    rho_c = parms$rho_c
    rho_l = parms$rho_l
    rho_o = parms$rho_o
    tau_parm = parms$tau_parm
    kappa_c = parms$kappa_c
    kappa_l = parms$kappa_l

    I = aggregates$I
    c_S = aggregates$c_S
    c_I = aggregates$c_I
    c_R = aggregates$c_R
    l_S = aggregates$l_S
    l_I = aggregates$l_I
    l_R = aggregates$l_R

    if(type=="S") { tau_j = rho_c*I*(c_S*c_I)^kappa_c + rho_l*I*(l_S*l_I)^kappa_l + rho_o*I }
    if(type=="I") { tau_j = rho_c*I*(c_I*c_I)^kappa_c + rho_l*I*(l_I*l_I)^kappa_l + rho_o*I }
    if(type=="R") { tau_j = rho_c*I*(c_R*c_I)^kappa_c + rho_l*I*(l_R*l_I)^kappa_l + rho_o*I }    
    # n_contacts = contacts(parms, aggregates)

    # P_I = tau_parm*n_contacts*I
    P_I = tau_parm*tau_j

    return(P_I)
}

contacts <- function(parms,choices) {
    rho_c = parms$rho_c
    rho_l = parms$rho_l
    rho_o = parms$rho_o
    kappa_c = parms$kappa_c
    kappa_l = parms$kappa_l

    c_S = choices$c_S
    c_I = choices$c_I
    l_S = choices$l_S
    l_I = choices$l_I

    n_contacts = rho_c*(c_S*c_I)^kappa_c + rho_l*(l_S*l_I)^kappa_l + rho_o

    return(n_contacts)
}

newly_infected <- function(parms,choices,sizes) {
    rho_c = parms$rho_c
    rho_l = parms$rho_l
    rho_o = parms$rho_o
    tau_parm = parms$tau_parm
    kappa_c = parms$kappa_c
    kappa_l = parms$kappa_l

    c_St = choices$c_S
    c_It = choices$c_I
    l_St = choices$l_S
    l_It = choices$l_I

    S_t = sizes$S
    I_t = sizes$I

    n_contacts = (rho_c*(c_St*c_It)^kappa_c + rho_l*(l_St*l_It)^kappa_l + rho_o)

    count = tau_parm*n_contacts*S_t*I_t

    # count = tau_parm*(rho_c*(c_St*S_t)*(c_It*I_t) + rho_l*(l_St*S_t)*(l_It*I_t) + rho_o*I_t*S_t)

    return(count)
}

# function to calculate the next SIR state
SIR_update <- function(parms,choices,aggregates) {
    rho_c = parms$rho_c
    rho_l = parms$rho_l
    rho_o = parms$rho_o
    recovery_rate = parms$pi_r
    case_fatality_rate = parms$pi_d
    tau_parm = parms$tau_parm
    kappa_c = parms$kappa_c
    kappa_l = parms$kappa_l

    c_S = choices$c_S
    c_I = choices$c_I
    l_S = choices$l_S
    l_I = choices$l_I

    S = aggregates$S
    I = aggregates$I
    R = aggregates$R
    D = aggregates$D

    parms_to_N_I = data.frame(rho_c = rho_c, rho_l = rho_l, rho_o = rho_o, tau_parm = tau_parm, kappa_c = kappa_c, kappa_l = kappa_l)
    choices_to_N_I = data.frame(c_S = c_S, c_I = c_I, l_S = l_S, l_I = l_I)
    SI_to_N_I = data.frame(S = S, I = I)

    T = max(newly_infected(parms_to_N_I, choices_to_N_I, SI_to_N_I),0)
    if(is.na(T)) { T = 0 }
    S_ = max(S - T,0)
    I_ = max((1 - recovery_rate - case_fatality_rate)*I + T,0)
    R_ = max(R + recovery_rate*I,0)
    D_ = max(D + case_fatality_rate*I,0)

    updated_SIR_state = data.frame(S=S_, I=I_, R=R_, D=D_)

    return(updated_SIR_state)
}

# function to generate important epi and econ series from solved policy functions and initial conditions
SIR_series <- function(parms,choices_list,aggregates,method="multilinear",grid_dfrm=NULL) {
    w = parms$w
    phi = parms$phi
    final.time = parms$final.time
    pi_r = parms$pi_r
    pi_d = parms$pi_d
    rho_c = parms$rho_c
    rho_l = parms$rho_l
    rho_o = parms$rho_o
    kappa_c = parms$kappa_c
    kappa_l = parms$kappa_l

    S_0 = aggregates$S
    I_0 = aggregates$I
    R_0 = aggregates$R
    D_0 = 1 - (S_0 + I_0 + R_0)

    state = data.frame(S=S_0,I=I_0,R=R_0)
    types = c("S","I","R")
    labor_supply = data.frame(S=NA,I=NA,R=NA)
    lifetime_utility = data.frame(S=NA,I=NA,R=NA)

    l_S0 = project_choices("S", state, choices_list$grid_list, choices_list$choices)
    l_I0 = project_choices("I", state, choices_list$grid_list, choices_list$choices)
    l_R0 = project_choices("R", state, choices_list$grid_list, choices_list$choices)

    output_series = data.frame(S = rep(S_0,length.out=final.time),
                                I = rep(I_0,length.out=final.time),
                                R = rep(R_0,length.out=final.time),
                                D = rep(D_0,length.out=final.time),
                                labor_S = rep(l_S0,length.out=final.time),
                                labor_I = rep(l_I0,length.out=final.time),
                                labor_R = rep(l_R0,length.out=final.time),
                                consumption_S = rep(w*l_S0,length.out=final.time),
                                consumption_I = rep(w*phi*l_I0,length.out=final.time),
                                consumption_R = rep(w*l_R0,length.out=final.time),
                                aggregate_labor_supply = rep((l_S0*S_0 + l_I0*I_0 + l_R0*R_0),length.out=final.time),
                                aggregate_consumption = rep((w*l_S0*S_0 + w*phi*l_I0*I_0 + w*l_R0*R_0),length.out=final.time)) 

    for(i in 2:final.time) {
        if(method=="multilinear") {
            for(j in 1:length(types)) {
                choices = choices_list$choices
                agent_type = types[j]
                labor_supply[j] = project_choices(agent_type, state, choices_list$grid_list, choices)
            }
        }

        l_S = labor_supply$S
        l_I = labor_supply$I

        c_S = w*l_S
        c_I = w*phi*l_I

        choices = data.frame(c_S = c_S, c_I = c_I, l_S = l_S, l_I = l_I)

        state = SIR_update(parms,choices,state)
        
        choice_vector = data.frame(
            labor_S = labor_supply$S, 
            labor_I = labor_supply$I, 
            labor_R = labor_supply$R,
            consumption_S = w*labor_supply$S,
            consumption_I = w*phi*labor_supply$I,
            consumption_R = w*labor_supply$R)
        aggregate_labor_supply = labor_supply$S*state$S + labor_supply$I*state$I + labor_supply$R*state$R
        aggregate_consumption = choice_vector$consumption_S*state$S + choice_vector$consumption_I*state$I + choice_vector$consumption_R*state$R

        output_series[i,] = c(state, choice_vector, aggregate_labor_supply, aggregate_consumption)
    }

    Reff_choices = data.frame(l_S=output_series$labor_S, l_I=output_series$labor_I, l_R=output_series$labor_R, c_S=output_series$consumption_S, c_I=output_series$consumption_I, c_R=output_series$consumption_R)
    Reff_SIR = data.frame(S=output_series$S, I=output_series$I, R=output_series$R, D=output_series$D)
    Reff = compute_R0(parms, Reff_choices, Reff_SIR)

    output_series = cbind(output_series, 
                    weekly_new_cases = newly_infected(parms, Reff_choices, Reff_SIR),
                    newly_recovered = pi_r*output_series$I,
                    newly_dead = pi_d*output_series$I,
                    newly_infected_per_I = newly_infected(parms, Reff_choices, Reff_SIR)/output_series$I,
                    prob_infection_S = newly_infected(parms, Reff_choices, Reff_SIR)/output_series$S,
                    infection_growth_rate = newly_infected(parms, Reff_choices, Reff_SIR)/output_series$I - pi_r - pi_d,
                    Reff = Reff,
                    consumption_contacts = rho_c*(output_series$consumption_S*output_series$consumption_I)^kappa_c,
                    labor_contacts = rho_l*(output_series$labor_S*output_series$labor_I)^kappa_l,
                    other_contacts = rho_o)
    output_series = cbind(output_series, total_cases = cumsum(output_series$weekly_new_cases), aggregate_hours_deviation = 100*(output_series$aggregate_labor_supply/output_series$aggregate_labor_supply[1] - 1) , aggregate_consumption_deviation = 100*(output_series$aggregate_consumption/output_series$aggregate_consumption[1] - 1) )

    return(output_series)
}


SIR_series_policy <- function(parms,choices_list,aggregates,policy_list) {
    w = parms$w
    phi = parms$phi
    final.time = parms$final.time
    pi_r = parms$pi_r
    pi_d = parms$pi_d
    rho_c = parms$rho_c
    rho_l = parms$rho_l
    rho_o = parms$rho_o
    kappa_c = parms$kappa_c
    kappa_l = parms$kappa_l

    S_0 = aggregates$S
    I_0 = aggregates$I
    R_0 = aggregates$R
    D_0 = 1 - (S_0 + I_0 + R_0)

    policy_switch = policy_list$switch
    policy_start = policy_list$start_date
    policy_stop = policy_list$end_date
    policy_labor_restriction = policy_list$labor_supply_level
    policy_state_dependent = policy_list$state_dependent

    state = data.frame(S=S_0,I=I_0,R=R_0)
    types = c("S","I","R")
    labor_supply = data.frame(S=NA,I=NA,R=NA)
    lifetime_utility = data.frame(S=NA,I=NA,R=NA)

    l_S0 = project_choices("S", state, choices_list$grid_list, choices_list$choices)
    l_I0 = project_choices("I", state, choices_list$grid_list, choices_list$choices)
    l_R0 = project_choices("R", state, choices_list$grid_list, choices_list$choices)

    restricted_level = l_S0*policy_labor_restriction

    output_series = data.frame(S = rep(S_0,length.out=final.time),
                                I = rep(I_0,length.out=final.time),
                                R = rep(R_0,length.out=final.time),
                                D = rep(D_0,length.out=final.time),
                                labor_S = rep(l_S0,length.out=final.time),
                                labor_I = rep(l_I0,length.out=final.time),
                                labor_R = rep(l_R0,length.out=final.time),
                                consumption_S = rep(w*l_S0,length.out=final.time),
                                consumption_I = rep(w*phi*l_I0,length.out=final.time),
                                consumption_R = rep(w*l_R0,length.out=final.time),
                                aggregate_labor_supply = rep((l_S0*S_0 + l_I0*I_0 + l_R0*R_0),length.out=final.time),
                                aggregate_consumption = rep((w*l_S0*S_0 + w*phi*l_I0*I_0 + w*l_R0*R_0),length.out=final.time)) 

    for(i in 2:final.time) {
        for(j in 1:length(types)) {
            choices = choices_list$choices
            agent_type = types[j]
            labor_supply[j] = project_choices(agent_type, state, choices_list$grid_list, choices)
            if((policy_switch==1)&(i>=policy_start)&(i<=policy_stop)) {
                labor_supply[j] = min(labor_supply[j],restricted_level)
            }
        }

        l_S = labor_supply$S
        l_I = labor_supply$I

        c_S = w*l_S
        c_I = w*phi*l_I

        choices = data.frame(c_S = c_S, c_I = c_I, l_S = l_S, l_I = l_I)

        state = SIR_update(parms,choices,state)
        choice_vector = data.frame(
            labor_S = labor_supply$S, 
            labor_I = labor_supply$I, 
            labor_R = labor_supply$R,
            consumption_S = w*labor_supply$S,
            consumption_I = w*phi*labor_supply$I,
            consumption_R = w*labor_supply$R)
        aggregate_labor_supply = labor_supply$S*state$S + labor_supply$I*state$I + labor_supply$R*state$R
        aggregate_consumption = choice_vector$consumption_S*state$S + choice_vector$consumption_I*state$I + choice_vector$consumption_R*state$R

        output_series[i,] = c(state, choice_vector, aggregate_labor_supply, aggregate_consumption)
    }

    Reff_choices = data.frame(l_S=output_series$labor_S, l_I=output_series$labor_I, l_R=output_series$labor_R, c_S=output_series$consumption_S, c_I=output_series$consumption_I, c_R=output_series$consumption_R)
    Reff_SIR = data.frame(S=output_series$S, I=output_series$I, R=output_series$R, D=output_series$D)
    Reff = compute_R0(parms, Reff_choices, Reff_SIR)

    output_series = cbind(output_series, 
                weekly_new_cases = newly_infected(parms, Reff_choices, Reff_SIR),
                newly_recovered = pi_r*output_series$I,
                newly_dead = pi_d*output_series$I,
                newly_infected_per_I = newly_infected(parms, Reff_choices, Reff_SIR)/output_series$I,
                prob_infection_S = newly_infected(parms, Reff_choices, Reff_SIR)/output_series$S,
                infection_growth_rate = newly_infected(parms, Reff_choices, Reff_SIR)/output_series$I - pi_r - pi_d,
                Reff = Reff,
                consumption_contacts = rho_c*(output_series$consumption_S*output_series$consumption_I)^kappa_c,
                labor_contacts = rho_l*(output_series$labor_S*output_series$labor_I)^kappa_l,
                other_contacts = rho_o)
    
    output_series = cbind(output_series, total_cases = cumsum(output_series$weekly_new_cases), aggregate_hours_deviation = 100*(output_series$aggregate_labor_supply/output_series$aggregate_labor_supply[1] - 1), aggregate_consumption_deviation = 100*(output_series$aggregate_consumption/output_series$aggregate_consumption[1] - 1))

    return(output_series)
}



# function to generate important epi and econ series from solved policy functions and initial conditions for impulse response plots
SIR_series_IRF <- function(parms,choices_list,aggregates,method="multilinear",grid_dfrm=NULL,impulse_parm=NULL,impulse_init=4,impulse_new=2,impulse_time=NULL, solved_values=NULL) {
    w = parms$w
    phi = parms$phi
    final.time = parms$final.time
    pi_r = parms$pi_r
    pi_d = parms$pi_d
    rho_c = impulse_init/(parms$Cbm)^2
    rho_l = parms$rho_l
    rho_o = parms$rho_o

    # S_0 = aggregates$S
    # I_0 = aggregates$I
    # R_0 = aggregates$R
    S_0 = aggregates[[1]]
    I_0 = aggregates[[2]]
    R_0 = aggregates[[3]]
    D_0 = 1 - (S_0 + I_0 + R_0)

    # state = c(S=S_0,I=I_0,R=R_0,benchmark_cons_contacts=impulse_init)
    state = c(S=S_0,I=I_0,R=R_0)
    types = c("S","I","R")
    labor_supply = data.frame(S=NA,I=NA,R=NA)
    lifetime_utility = data.frame(S=NA,I=NA,R=NA)

    l_S = solved_values$labor_supply_S[which(solved_values$benchmark_cons_contacts==impulse_init)]
    l_I = solved_values$labor_supply_I[which(solved_values$benchmark_cons_contacts==impulse_init)]
    l_R = solved_values$labor_supply_R[which(solved_values$benchmark_cons_contacts==impulse_init)]

    l_S_fn = ipol(l_S, grid=grid_list, method="multilinear")
    l_I_fn = ipol(l_I, grid=grid_list, method="multilinear")
    l_R_fn = ipol(l_R, grid=grid_list, method="multilinear")
    l_S0 = l_S_fn(state)
    l_I0 = l_I_fn(state)
    l_R0 = l_R_fn(state)

    output_series = data.frame(S = rep(S_0,length.out=final.time),
                                I = rep(I_0,length.out=final.time),
                                R = rep(R_0,length.out=final.time),
                                D = rep(D_0,length.out=final.time),
                                labor_S = rep(l_S0,length.out=final.time),
                                labor_I = rep(l_I0,length.out=final.time),
                                labor_R = rep(l_R0,length.out=final.time),
                                consumption_S = rep(w*l_S0,length.out=final.time),
                                consumption_I = rep(w*phi*l_I0,length.out=final.time),
                                consumption_R = rep(w*l_R0,length.out=final.time),
                                aggregate_labor_supply = rep((l_S0*S_0 + l_I0*I_0 + l_R0*R_0),length.out=final.time),
                                aggregate_consumption = rep((w*l_S0*S_0 + w*phi*l_I0*I_0 + w*l_R0*R_0),length.out=final.time)) 

    print(l_S_fn(c(0.99,0.01,0)))
    print(parms)

    choices = choices_list$choices
    for(i in 2:(impulse_time-1)) {
        if(method=="multilinear") {
            for(j in 1:length(types)) {
                # if(i>=impulse_time) state[4] = impulse_new
                
                agent_type = types[j]
                # print(state)
                state = as.numeric(state)
                if(agent_type=="S") labor_supply[j] = l_S_fn(state)
                if(agent_type=="I") labor_supply[j] = l_I_fn(state)
                if(agent_type=="R") labor_supply[j] = l_R_fn(state)
                # print(paste0("worked ",i))
                # print(labor_supply$S)
            }
            # print("Collect 200")
        }

        l_S = labor_supply$S
        l_I = labor_supply$I

        # print("Collect 400")

        c_S = w*l_S
        c_I = w*phi*l_I

        choices = data.frame(c_S = c_S, c_I = c_I, l_S = l_S, l_I = l_I)

        # state_old = data.frame(S = state[1], I = state[2], R = state[3], benchmark_cons_contacts = state[4])
        # state = SIR_update(parms,choices,state_old)
        # print(state)
        state = data.frame(S = state[1], I = state[2], R = state[3])
        state = SIR_update(parms,choices,state)
        state = state[1:3]
        # state[4] = state_old[4]
        # print(state)
        choice_vector = data.frame(
            labor_S = labor_supply$S, 
            labor_I = labor_supply$I, 
            labor_R = labor_supply$R,
            consumption_S = w*labor_supply$S,
            consumption_I = w*phi*labor_supply$I,
            consumption_R = w*labor_supply$R)
        aggregate_labor_supply = labor_supply$S*state$S + labor_supply$I*state$I + labor_supply$R*state$R
        aggregate_consumption = choice_vector$consumption_S*state$S + choice_vector$consumption_I*state$I + choice_vector$consumption_R*state$R

        output_series[i,] = c(state, choice_vector, aggregate_labor_supply, aggregate_consumption)
    }

    l_S = solved_values$labor_supply_S[which(solved_values$benchmark_cons_contacts==impulse_new)]
    l_I = solved_values$labor_supply_I[which(solved_values$benchmark_cons_contacts==impulse_new)]
    l_R = solved_values$labor_supply_R[which(solved_values$benchmark_cons_contacts==impulse_new)]
    l_S_fn = ipol(l_S, grid=grid_list, method="multilinear")
    l_I_fn = ipol(l_I, grid=grid_list, method="multilinear")
    l_R_fn = ipol(l_R, grid=grid_list, method="multilinear")
    parms$rho_c = impulse_new/(parms$Cbm)^2

    print(l_S_fn(c(0.99,0.01,0)))

    for(i in impulse_time:final.time) {
        if(method=="multilinear") {
            for(j in 1:length(types)) {
                
                agent_type = types[j]
                state = as.numeric(state)
                if(agent_type=="S") labor_supply[j] = l_S_fn(state)
                if(agent_type=="I") labor_supply[j] = l_I_fn(state)
                if(agent_type=="R") labor_supply[j] = l_R_fn(state)
            }
        }

        l_S = labor_supply$S
        l_I = labor_supply$I

        c_S = w*l_S
        c_I = w*phi*l_I

        choices = data.frame(c_S = c_S, c_I = c_I, l_S = l_S, l_I = l_I)

        state = data.frame(S = state[1], I = state[2], R = state[3])
        state = SIR_update(parms,choices,state)
        state = state[1:3]
        choice_vector = data.frame(
            labor_S = labor_supply$S, 
            labor_I = labor_supply$I, 
            labor_R = labor_supply$R,
            consumption_S = w*labor_supply$S,
            consumption_I = w*phi*labor_supply$I,
            consumption_R = w*labor_supply$R)
        aggregate_labor_supply = labor_supply$S*state$S + labor_supply$I*state$I + labor_supply$R*state$R
        aggregate_consumption = choice_vector$consumption_S*state$S + choice_vector$consumption_I*state$I + choice_vector$consumption_R*state$R

        output_series[i,] = c(state, choice_vector, aggregate_labor_supply, aggregate_consumption)
    }

    Reff_choices = data.frame(l_S=output_series$labor_S, l_I=output_series$labor_I, l_R=output_series$labor_R, c_S=output_series$consumption_S, c_I=output_series$consumption_I, c_R=output_series$consumption_R)
    Reff_SIR = data.frame(S=output_series$S, I=output_series$I, R=output_series$R, D=output_series$D)
    Reff = compute_R0(parms, Reff_choices, Reff_SIR)

    output_series = cbind(output_series, 
                    weekly_new_cases = newly_infected(parms, Reff_choices, Reff_SIR),
                    newly_recovered = pi_r*output_series$I,
                    newly_dead = pi_d*output_series$I,
                    newly_infected_per_I = newly_infected(parms, Reff_choices, Reff_SIR)/output_series$I,
                    prob_infection_S = newly_infected(parms, Reff_choices, Reff_SIR)/output_series$S,
                    infection_growth_rate = newly_infected(parms, Reff_choices, Reff_SIR)/output_series$I - pi_r - pi_d,
                    Reff = Reff,
                    consumption_contacts = rho_c*output_series$consumption_S*output_series$consumption_I,
                    labor_contacts = rho_l*output_series$labor_S*output_series$labor_I,
                    other_contacts = rho_o)
    output_series = cbind(output_series, total_cases = cumsum(output_series$weekly_new_cases), aggregate_hours_deviation = 100*(output_series$aggregate_labor_supply/output_series$aggregate_labor_supply[1] - 1) , aggregate_consumption_deviation = 100*(output_series$aggregate_consumption/output_series$aggregate_consumption[1] - 1) )

    return(output_series)
}

# steady-state utility function
ss_utility <- function(args, discount_factor) {
    c = args$c
    l = args$l
    C = args$C
    L = args$L
    gamma_c = args$gamma_c
    gamma_l = args$gamma_l
    rho_c = args$rho_c
    rho_l = args$rho_l
    s = args$s
    risk_aversion = args$risk_aversion
    alpha_U = args$alpha_U
    r = (1-discount_factor)/discount_factor
    Lbar = args$Lbar

    num = total_utility(c,l,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha_U,Lbar)
    den = 1-discount_factor

    # value = num + discount_factor*(num/den)
    value = num/den

    return(value)
}

composite_good <- function(c,l,C,L,gamma_c,gamma_l,rho_c,rho_l,s,alpha,Lbar=24) {
    consumption_piece = (alpha^(1/s))*c^((s-1)/s)
    leisure_piece = ((1-alpha)^(1/s))*(Lbar-l)^((s-1)/s)
    private_piece = (consumption_piece + leisure_piece)^(s/(s-1))
    contact_piece = (1 + gamma_c*rho_c*(c*C)^0.25 + gamma_l*rho_l*(l*L)^0.25)
    composite_good = private_piece*contact_piece

    return(composite_good)
}

total_utility <- function(c,l,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha_U,Lbar=24) {
    
    good = composite_good(c,l,C,L,gamma_c,gamma_l,rho_c,rho_l,s,alpha_U,Lbar)

    if(risk_aversion==1) { value =  log(good) }
    if(risk_aversion!=1) { value = (good^(1-risk_aversion) - 1)/(1 - risk_aversion)}

    return(value)
}

# wrapper over static objective function to use in calibration
TU_wrapper <- function(choice,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha,w,p=1,Lbar=24,nonlabor_income=0,phi=1) {
    l = choice
    c = (phi*w*choice + nonlabor_income)/p
    L = l
    C = c
    value = total_utility(c,l,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha,Lbar)
    return(value)
}

# wrapper over planner's objective function to use in calibration
planner_TU_wrapper <- function(choices,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha,w,p=1,Lbar=24,nonlabor_income=0,phi,state=c(1-1e-06,1e-06,0)) {
    l_S = choices[1]
    c_S = (w*choice[1] + nonlabor_income)/p
    l_I = choices[2]
    c_I = (w*phi*choice[2] + nonlabor_income)/p
    l_R = choices[3]
    c_R = (w*choice[3] + nonlabor_income)/p
    
    S = state[1]
    I = state[2]
    R = state[3]

    L = S*l_S + I*l_I + R*l_R
    C = S*c_S + I*c_I + R*c_R
    value = total_utility(c,l,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha,Lbar)
    return(value)
}

## wrapper over dynamic objective function to verify calibration works & use in calibration
dynamic_TU_wrapper <- function(choice, parms, benchmarks, state) {
    w = benchmarks$w
    p = benchmarks$p
    rho_c = benchmarks$rho_c
    rho_l = benchmarks$rho_l
    risk_aversion = benchmarks$risk_aversion
    discount_factor = benchmarks$discount_factor
    Lbm = benchmarks$Lbm
    Cbm = benchmarks$Cbm

    s = parms[1]
    alpha = parms[2]
    gamma_c = parms[3]
    gamma_l = parms[4]

    S = state$S
    I = state$I
    R = state$R

    l = choice
    c = w*choice/p
    L = S*l + I*0.8*l + R*Lbm
    C = S*c + I*0.8*c + R*Cbm

    aggregates = data.frame(I=state$I, c_S=c, l_S=l, c_I=0.8*c, l_I=0.8*l, c_R=Cbm, l_R=Lbm)

    P_I = tau_t(type="S", benchmarks, aggregates)

    cv_S = total_utility(c,l,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha,Lbar)/(1-discount_factor)
    cv_I = total_utility(0.8*c,0.8*l,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha,Lbar)/(1-discount_factor)

    # if(type=="S") {
        value = total_utility(c,l,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha,Lbar) + P_I*cv_I + (1-P_I)*cv_S
    # }

    return(value)
}

project_SIR_outcome <- function(SIR_state, grid_list, target) {
    state = as.matrix(SIR_state)
    projection = ipol(target, grid=grid_list, method="multilinear")
    value = projection(state)

    return(value)
}

# continuation value interpolation created from a grid of function values and evaluation nodes
continuation_value <- function(updated_SIR_state, grid_list, contval_dfrm) {
    state = c(updated_SIR_state$S, updated_SIR_state$I, updated_SIR_state$R)
    cv_fn = ipol(contval_dfrm, grid=grid_list, method="multilinear")
    value = cv_fn(state)

    return(value)
}

# FOC for the measure-zero S-type
m0_S <- function(choice, exog_parms, others_choices, SIR_state, grid_list, contval_list) {
   # structural utility parameters
    gamma_c = exog_parms$gamma_c
    gamma_l = exog_parms$gamma_l
    risk_aversion = exog_parms$risk_aversion
    s = exog_parms$s
    alpha_U = exog_parms$alpha_U
    Lbar = exog_parms$Lbar
    # structural contact parameters
    rho_c = exog_parms$rho_c
    rho_l = exog_parms$rho_l
    rho_o = exog_parms$rho_o
    tau_parm = exog_parms$tau_parm
    pi_r = exog_parms$pi_r
    pi_d = exog_parms$pi_d
    kappa_c = exog_parms$kappa_c
    kappa_l = exog_parms$kappa_l
    # other structural parameters
    discount_factor = exog_parms$discount_factor
    phi = exog_parms$phi
    w = exog_parms$w
    Lbm = exog_parms$Lbm
    # m0type = exog_parms$m0type

    # current SIR state
    S = SIR_state$S
    I = SIR_state$I
    R = SIR_state$R

    # choices
    l = choice[[1]]
    l_I = as.numeric(others_choices$I)
    l_R = as.numeric(others_choices$R)
    c = w*l
    c_I = w*phi*l_I
    c_R = w*l_R

    choices_dfrm = data.frame(c_S = c, c_I = c_I, l_S = l, l_I = l_I) 
    SIR_update_parms = data.frame(rho_c, rho_l, rho_o, pi_r, pi_d, phi, tau_parm, kappa_c, kappa_l)

    SIR_ = SIR_update(parms=SIR_update_parms, choices=choices_dfrm, aggregates=SIR_state)

    C = c*S + c_I*I + c_R*R
    L = l*S + l_I*I + l_R*R

    # CES aggregation of consumption and labor
    Cbar = composite_good(c,l,C,L,gamma_c,gamma_l,rho_c,rho_l,s,alpha_U,Lbar)
    contact_piece = (1 + gamma_c*rho_c*(c*C)^0.25 + gamma_l*rho_l*(l*L)^0.25)

    # MU of consumption
    #u_c = (Cbar^(-1-risk_aversion))*(alpha_U/c)^(1/s)
    dCbar_dc = (Cbar^(-1))*((alpha_U/c)^(1/s))*contact_piece + (Cbar*0.25*gamma_c*rho_c*(c*C)^(0.25))/c
    u_c = (Cbar^(-risk_aversion))*dCbar_dc

    # MU of labor
    #u_l = -(Cbar^(-1-risk_aversion))*((1-alpha_U)/(Lbar-l))^(1/s)
    dCbar_dl = -(Cbar^(-1))*(((1-alpha_U)/(Lbar-l))^(1/s))*contact_piece + (Cbar*0.25*gamma_l*rho_l*(l*L)^(0.25))/l
    u_l = (Cbar^(-risk_aversion))*dCbar_dl

    # derivative of infection risk wrt S-type consumption
    # P_I_c = tau_parm*rho_c*c_I*I
    P_I_c = tau_parm*rho_c*I*(c_I^kappa_c)*kappa_c*c^(kappa_c-1)

    # derivative of infection risk wrt S-type labor
    # P_I_l = tau_parm*rho_l*l_I*I
    P_I_l = tau_parm*rho_l*I*(l_I^kappa_l)*kappa_l*l^(kappa_l-1)

    # probability gap for discrete-choice measure 0
    # prob_gap = tau_parm*Lbm*(w*rho_c*c_I + rho_l*l_I)*I

    # U^S
    U_S = continuation_value(SIR_, grid_list, contval_list$S)

    # U^I
    U_I = continuation_value(SIR_, grid_list, contval_list$I)

    # FOC
    FOC = u_c + u_l/w - discount_factor*(P_I_c + P_I_l/w)*(U_S - U_I)

    return(FOC)
}

m0_component_plot <- function(exog_parms, SIR_state, grid_list, contval_list, type="FOC") {
   # structural utility parameters
    gamma_c = exog_parms$gamma_c
    gamma_l = exog_parms$gamma_l
    risk_aversion = exog_parms$risk_aversion
    s = exog_parms$s
    alpha_U = exog_parms$alpha_U
    # structural contact parameters
    rho_c = exog_parms$rho_c
    rho_l = exog_parms$rho_l
    rho_o = exog_parms$rho_o
    tau_parm = exog_parms$tau_parm
    pi_r = exog_parms$pi_r
    pi_d = exog_parms$pi_d
    # other structural parameters
    discount_factor = exog_parms$discount_factor
    phi = exog_parms$phi
    max_labor_supply = exog_parms$max_labor_supply
    w = exog_parms$w
    y = exog_parms$nonlabor_income
    magic_time = exog_parms$magic_time_income
    Lbar = exog_parms$Lbar

   # SIR_state = data.frame(S=1,I=0,R=0)

    # current SIR state
    S = SIR_state$S
    I = SIR_state$I
    R = SIR_state$R

    l = seq(from=0, to=8.1, length.out=100)
    l_I = 7.8
    l_R = 8
    c = w*l + y
    c_I = w*phi*l_I + y
    c_R = w*l_R + y

    choices_dfrm = data.frame(c_S = c, c_I = c_I, l_S = l, l_I = l_I) 
    SIR_update_parms = data.frame(rho_c, rho_l, rho_o, pi_r, pi_d, phi, tau_parm)

    SIR_ = data.frame(matrix(-1,nrow=length(l),ncol=4))
    colnames(SIR_) = c("S","I","R","D")
    Cbar = rep(0,length=length(l))
    u_c = rep(0,length=length(l))
    u_l = rep(0,length=length(l))
    U_S = rep(0,length=length(l))
    U_I = rep(0,length=length(l))
    FOC = rep(0,length=length(l))
    TU = rep(0,length=length(l))
    percent_change_U = rep(0,length=length(l))
    n_contacts = rep(0,length=length(l))

    for(i in 1:length(l)) {

        SIR_[i,] = SIR_update(parms=SIR_update_parms, choices=choices_dfrm[i,], aggregates=SIR_state)

        C = c*S + c_I*I + c_R*R
        L = l*S + l_I*I + l_R*R

        # CES aggregation of consumption and labor
        Cbar[i] = composite_good(c[i],l[i],C[i],L[i],gamma_c,gamma_l,rho_c,rho_l,s,alpha_U,Lbar)
        # Cbar_1[i] = composite_good(w*8,8,C[i],L[i],gamma_c,gamma_l,rho_c,rho_l,s,alpha_U)
        # Cbar_0[i] = composite_good(0,0,C[i],L[i],gamma_c,gamma_l,rho_c,rho_l,s,alpha_U)

        # MU of consumption
        u_c[i] = (Cbar[i]^(-1-risk_aversion))*(alpha_U/c[i])^(1/s)

        # MU of labor
        u_l[i] = -(Cbar[i]^(-1-risk_aversion))*((1-alpha_U)/(Lbar-l[i] + magic_time))^(1/s)

        # TU of consumption as worker
        # u_1 = (Cbar_1^(1 - risk_aversion) - 1)/(1-risk_aversion)

        # # # TU of consumption as worker
        # u_0 = (Cbar_0^(1 - risk_aversion) - 1)/(1-risk_aversion)

        # derivative of infection risk wrt S-type consumption
        P_I_c = tau_parm*rho_c*c_I*I

        # derivative of infection risk wrt S-type labor
        P_I_l = tau_parm*rho_l*l_I*I

        # number of contacts
        n_contacts[i] = rho_c*c[i]*c_I + rho_l*l[i]*l_I + rho_o

        # difference in infection probabilities under discrete labor supply choice
        #prob_gap = tau_parm*8*(w*rho_c*c_I + rho_l*l_I)*I

        # U^S
        U_S[i] = continuation_value(as.data.frame(SIR_[i,]), grid_list, contval_list$lifetime_utility_S)

        # U^I
        U_I[i] = continuation_value(SIR_[i,], grid_list, contval_list$lifetime_utility_I)

        # FOC
        FOC[i] = u_c[i] + u_l[i]/w - discount_factor*(P_I_c + P_I_l/w)*(U_S[i]- U_I[i])
        # FOC_2[i] = u_1[i] - u_0[i] - discount_factor*prob_gap*(U_S[i] - U_I[i])

        TU[i] = total_utility(c[i],l[i],C[i],L[i],gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha_U,Lbar)

        percent_change_U[i] = (u_c[i] + u_l[i]/w)/TU[i]

    }

    I_grid = seq(from=0,to=0.2,length.out=100)
    U_S_over_I = rep(0,length=length(I_grid)) 
    U_I_over_I = rep(0,length=length(I_grid)) 
    n_contacts_over_I = rep(0,length=length(I_grid)) 
    P_I_c_over_I = rep(0,length=length(I_grid)) 
    P_I_l_over_I = rep(0,length=length(I_grid)) 
    prob_sum_over_I = rep(0,length=length(I_grid)) 

    SIR_over_I = data.frame(matrix(-1,nrow=length(l),ncol=4))
    colnames(SIR_over_I) = c("S","I","R","D")
    choices_dfrm_small = data.frame(S = 0, I = l_I, R = l_R) 
    choices_dfrm = data.frame(c_S = U_S_over_I, c_I = c_I, l_S = U_S_over_I, l_I = l_I) 
    for(i in 1:length(I_grid)) {

        SIR_state$I = I_grid[i]
        SIR_state$R = 0
        SIR_state$S = 1 - SIR_state$I
        result = tryCatch(uniroot(m0_S, lower=0, upper=Lbar, exog_parms = exog_parms, others_choices = choices_dfrm_small, SIR_state = SIR_state, grid_list = grid_list, contval_list = contval_list), error = function(err) {
                object = data.frame(root=0)
            })

        choices_dfrm$l_S[i] = result$root
        choices_dfrm$c_S[i] = result$root*w + y

        # SIR state update over I
        SIR_over_I[i,] = SIR_update(parms=SIR_update_parms, choices=choices_dfrm[i,], aggregates=SIR_state)
        
        # RHS over I
        U_S_over_I[i] = continuation_value(SIR_over_I[i,], grid_list, contval_list$lifetime_utility_S)
        U_I_over_I[i] = continuation_value(SIR_over_I[i,], grid_list, contval_list$lifetime_utility_I)
        P_I_c_over_I[i] = tau_parm*rho_c*c_I*I_grid[i]
        P_I_l_over_I[i] = tau_parm*rho_l*l_I*I_grid[i]

        prob_sum_over_I[i] = (P_I_c_over_I[i] + P_I_l_over_I[i]/w)*I_grid[i]

        # number of contacts over I
        n_contacts_over_I[i] = rho_c*choices_dfrm$c_S[i]*c_I + rho_l*choices_dfrm$l_S[i]*l_I + rho_o
    }

#    print(percent_change_U)

    ls_8 = which.min( (l-8)^2 )
    ls_6 = which.min( (l-6)^2 )
    ls_4 = which.min( (l-4)^2 )
    print(paste0("Going from 8 to 6: ", (TU[ls_8] - TU[ls_6])/TU[ls_8]))
    print(paste0("Going from 6 to 4: ", (TU[ls_6] - TU[ls_4])/TU[ls_6]))
    print(paste0("Going from 8 to 4: ", (TU[ls_8] - TU[ls_4])/TU[ls_8]))

    lines = data.frame(l=l, u_c = u_c, u_l = u_l, U_S = U_S, U_I = U_I, LHS = u_c + u_l/w, RHS = discount_factor*(P_I_c + P_I_l/w)*(U_S- U_I), TU = TU, pcU = percent_change_U, S_prime = SIR_$S, I_prime = SIR_$I, U_S = U_S, contacts = n_contacts)
    lines_over_I = data.frame(I = I_grid, S_prime = SIR_over_I$S, I_prime = SIR_over_I$I, U_S = U_S_over_I, RHS_over_I = discount_factor*prob_sum_over_I*(U_S_over_I- U_I_over_I), contacts = n_contacts_over_I)

    if(type=="FOC"){
        plot = ggplot(data = lines, aes(x=l)) + geom_line(aes(y = LHS), size = 1) +
                                        geom_line(aes(y = RHS), size = 1, color="red") +
                                        geom_hline(aes(yintercept=0), size=0.5, color="darkgray") +
                                        ylim(c(0,0.0009)) +
                                        theme_bw() + ylab("LHS (black) and RHS (red) of FOC") + xlab("S-type labor supply") + 
                                        ggtitle(paste0("SIR = (",S,",",I,",",R,"), rho_c = ", rho_c,"\n alpha = ", round(alpha_U,3), ", sigma = ", round(s,3), ", Lbar = ", Lbar, " hours. Wage = $", round(w,3)
                                            )
                                            )
    }

    if(type=="MBMC"){
        plot = ggplot(data = lines, aes(x=max(contacts) - contacts) ) + geom_line(aes(y = LHS), size = 1) +
                                        geom_line(aes(y = RHS), size = 1, color="red") +
                                        geom_hline(aes(yintercept=0), size=0.5, color="darkgray") +
                                        ylim(c(0,0.0009)) +
                                        theme_bw() + ylab("PMC_S (black) and PMB_S (red) of contact reduction") + xlab("S-type's contacts reduced") + 
                                        ggtitle(paste0("SIR = (",S,",",I,",",R,"), RA = ", risk_aversion,"\n alpha = ", round(alpha_U,3), ", sigma = ", round(s,3), ", Lbar = ", Lbar, " hours. Wage = $", round(w,3)
                                            )
                                            )
    }

    if(type=="total_utility"){
        plot = ggplot(data = lines, aes(x=l))  + geom_line(aes(y = TU), size = 1) +
                                        theme_bw() + ylab("Total utility") + xlab("S-type labor supply") + 
                                        ggtitle(paste0("SIR = (",S,",",I,",",R,"), rho_c = ", rho_c,"\n Daily nonlabor income: $", y
                                            )
                                            )
    }

    if(type=="utility_perc_change"){
        plot = ggplot(data = lines, aes(x=l))  + geom_line(aes(y = pcU), size = 1) +
                                        theme_bw() + ylab("% change in present utility") + xlab("S-type labor supply") + 
                                        xlim(c(6,8)) +
                                        ggtitle(paste0("SIR = (",S,",",I,",",R,"), rho_c = ", rho_c,"\n alpha = ", round(alpha_U,3), ", sigma = ", round(s,3), ", Lbar = ", Lbar, " hours."
                                            )
                                            )
    }

    if(type=="S_prime_ls"){
        plot = ggplot(data = lines, aes(x=l)) + geom_line(aes(y = S_prime), size = 1) +
                                        theme_bw() + ylab("S_{t+1}") + xlab("S-type labor supply") + 
                                        ggtitle(paste0("SIR = (",S,",",I,",",R,"), rho_c = ", rho_c,"\n alpha = ", round(alpha_U,3), ", sigma = ", round(s,3), ", Lbar = ", Lbar, " hours."
                                            )
                                            )
    }

    if(type=="I_prime_ls"){
        plot = ggplot(data = lines, aes(x=l)) + geom_line(aes(y = I_prime), size = 1) +
                                        theme_bw() + ylab("I_{t+1}") + xlab("S-type labor supply") + 
                                        ggtitle(paste0("SIR = (",S,",",I,",",R,"), rho_c = ", rho_c,"\n alpha = ", round(alpha_U,3), ", sigma = ", round(s,3), ", Lbar = ", Lbar, " hours."
                                            )
                                            )
    }

    if(type=="U_S"){
        plot = ggplot(data = lines, aes(x=l)) + geom_line(aes(y = U_S), size = 1) +
                                        theme_bw() + ylab("U_S") + xlab("S-type labor supply") + 
                                        ggtitle(paste0("SIR = (",S,",",I,",",R,"), rho_c = ", rho_c,"\n alpha = ", round(alpha_U,3), ", sigma = ", round(s,3), ", Lbar = ", Lbar, " hours."
                                            )
                                            )
    }

    if(type=="contacts"){
        plot = ggplot(data = lines, aes(x=l)) + geom_line(aes(y = contacts), size = 1) +
                                        theme_bw() + ylab("Number of contacts") + xlab("S-type labor supply") + 
                                        ggtitle(paste0("SIR = (",S,",",I,",",R,"), rho_c = ", rho_c,"\n alpha = ", round(alpha_U,3), ", sigma = ", round(s,3), ", Lbar = ", Lbar, " hours."
                                            )
                                            )
    }

    if(type=="RHS"){
        plot = ggplot(data = lines, aes(x=l)) + geom_line(aes(y = RHS), size = 1) +
                                        theme_bw() + ylab("RHS of FOC") + xlab("S-type labor supply") + 
                                        ggtitle(paste0("SIR = (",S,",",I,",",R,"), rho_c = ", rho_c,"\n alpha = ", round(alpha_U,3), ", sigma = ", round(s,3), ", Lbar = ", Lbar, " hours."
                                            )
                                            )
    }

    if(type=="RHS_over_I"){
        plot = ggplot(data = lines_over_I, aes(x=I)) + geom_line(aes(y = RHS_over_I), size = 1) +
                                        theme_bw() + ylab("RHS of FOC") + xlab("I_t") + 
                                        ggtitle(paste0("rho_c = ", rho_c,"\n alpha = ", round(alpha_U,3), ", sigma = ", round(s,3), ", Lbar = ", Lbar, " hours."
                                            )
                                            )
    }

    if(type=="contacts_over_I"){
        plot = ggplot(data = lines_over_I, aes(x=I)) + geom_line(aes(y = contacts), size = 1) +
                                        theme_bw() + ylab("Number of contacts") + xlab("I_t") + 
                                        ggtitle(paste0("rho_c = ", rho_c,"\n alpha = ", round(alpha_U,3), ", sigma = ", round(s,3), ", Lbar = ", Lbar, " hours."
                                            )
                                            )
    }

    if(type=="S_prime_I"){
        plot = ggplot(data = lines_over_I, aes(x=I)) + geom_line(aes(y = S_prime), size = 1) +
                                        theme_bw() + ylab("S_{t+1}") + xlab("I_t") + 
                                        ggtitle(paste0("rho_c = ", rho_c,"\n alpha = ", round(alpha_U,3), ", sigma = ", round(s,3), ", Lbar = ", Lbar, " hours."
                                            )
                                            )
    }

    if(type=="I_prime_I"){
        plot = ggplot(data = lines_over_I, aes(x=I)) + geom_line(aes(y = I_prime), size = 1) +
                                        theme_bw() + ylab("I_{t+1}") + xlab("I_t") + 
                                        ggtitle(paste0("rho_c = ", rho_c,"\n alpha = ", round(alpha_U,3), ", sigma = ", round(s,3), ", Lbar = ", Lbar, " hours."
                                            )
                                            )
    }

    if(type=="U_S_over_I"){
        plot = ggplot(data = lines_over_I, aes(x=l)) + geom_line(aes(y = U_S), size = 1) +
                                        theme_bw() + ylab("U_S") + xlab("I_t") + 
                                        ggtitle(paste0("rho_c = ", rho_c,"\n alpha = ", round(alpha_U,3), ", sigma = ", round(s,3), ", Lbar = ", Lbar, " hours."
                                            )
                                            )
    }

    return(plot)
}

# pre-value function of type j. they are aware of the SIR model and are able to compute future states.
pre_value_fn <- function(choice, exog_parms, others_choices, SIR_state, agent_type, grid_list, contval_list) {
    # structural utility parameters
    gamma_c = exog_parms$gamma_c
    gamma_l = exog_parms$gamma_l
    risk_aversion = exog_parms$risk_aversion
    s = exog_parms$s
    alpha_U = exog_parms$alpha_U
    Lbar = exog_parms$Lbar
    kappa_c = exog_parms$kappa_c
    kappa_l = exog_parms$kappa_l
    # structural contact parameters
    rho_c = exog_parms$rho_c
    rho_l = exog_parms$rho_l
    rho_o = exog_parms$rho_o
    tau_parm = exog_parms$tau_parm
    pi_r = exog_parms$pi_r
    pi_d = exog_parms$pi_d
    # other structural parameters
    discount_factor = exog_parms$discount_factor
    phi = exog_parms$phi
    max_labor_supply = exog_parms$max_labor_supply
    w = exog_parms$w
    U_D = exog_parms$udeath

    # current SIR state
    S = SIR_state$S
    I = SIR_state$I
    R = SIR_state$R

    l = choice

    l_S = as.numeric(others_choices$S)[1]
    l_I = as.numeric(others_choices$I)[1]
    l_R = as.numeric(others_choices$R)[1]

    if(agent_type=="S") { 
        l_S = l
        c = w*l
    }
    if(agent_type=="I") { 
        l_I = l 
        c = w*phi*l
    }
    if(agent_type=="R") { 
        l_R = l
        c = w*l
    }

    c_S = w*l_S
    c_I = w*phi*l_I
    c_R = w*l_R
    C = c_S*S + c_I*I + c_R*R
    L = l_S*S + l_I*I + l_R*R

    choices_dfrm = data.frame(c_S = c_S, c_I = c_I, l_S = l_S, l_I = l_I) 
    SIR_update_parms = data.frame(rho_c, rho_l, rho_o, pi_r, pi_d, phi, tau_parm, kappa_c, kappa_l)

    SIR_ = SIR_update(parms=SIR_update_parms, choices=choices_dfrm, aggregates=SIR_state)

    parms = data.frame(rho_c=rho_c, rho_l=rho_l, rho_o=rho_o, tau_parm=tau_parm, kappa_c=kappa_c, kappa_l=kappa_l)
    aggregates = data.frame(I=I, c_S=c_S, c_I=c_I, c_R=c_R, l_S=l_S, l_I=l_I, l_R=l_R)
    tau_jt = tau_t(agent_type, parms, aggregates)

    # print(aggregates)
    # print(tau_jt)

    if(agent_type=="S") { 
    value = total_utility(c,l,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha_U,Lbar) + discount_factor*((1-tau_jt)*continuation_value(SIR_, grid_list, contval_list$S)+ tau_jt*continuation_value(SIR_, grid_list, contval_list$I))
    }
    if(agent_type=="I") { 
    value = total_utility(c,l,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha_U,Lbar) + discount_factor*((1-pi_r-pi_d)*continuation_value(SIR_, grid_list, contval_list$I) + pi_r*continuation_value(SIR_, grid_list, contval_list$R)) + pi_d*U_D
    }
    if(agent_type=="R") { 
    value = total_utility(c,l,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha_U,Lbar) + discount_factor*(continuation_value(SIR_, grid_list, contval_list$R))
    }

    return( value )
}

# pre-value function of the planner.
pre_value_fn_planner <- function(choices, exog_parms, SIR_state, grid_list, contval_list) {
    # structural utility parameters
    gamma_c = exog_parms$gamma_c
    gamma_l = exog_parms$gamma_l
    risk_aversion = exog_parms$risk_aversion
    s = exog_parms$s
    alpha_U = exog_parms$alpha_U
    Lbar = exog_parms$Lbar
    kappa_c = exog_parms$kappa_c
    kappa_l = exog_parms$kappa_l
    # structural contact parameters
    rho_c = exog_parms$rho_c
    rho_l = exog_parms$rho_l
    rho_o = exog_parms$rho_o
    tau_parm = exog_parms$tau_parm
    pi_r = exog_parms$pi_r
    pi_d = exog_parms$pi_d
    # other structural parameters
    discount_factor = exog_parms$discount_factor
    phi = exog_parms$phi
    Lbar = exog_parms$Lbar
    w = exog_parms$w
    U_D = exog_parms$udeath

    plannertype = exog_parms$plannertype

    # current SIR state
    S = SIR_state$S
    I = SIR_state$I
    R = SIR_state$R

    l_S = choices[[1]]
    l_I = choices[[2]]
    l_R = choices[[3]]

    c_S = w*l_S
    c_I = w*phi*l_I
    c_R = w*l_R
    C = c_S*S + c_I*I + c_R*R
    L = l_S*S + l_I*I + l_R*R

    choices_dfrm = data.frame(c_S = c_S, c_I = c_I, l_S = l_S, l_I = l_I) 
    SIR_update_parms = data.frame(rho_c, rho_l, rho_o, pi_r, pi_d, phi, tau_parm, kappa_c, kappa_l)

    SIR_ = SIR_update(parms=SIR_update_parms, choices=choices_dfrm, aggregates=SIR_state)

    parms = data.frame(rho_c=rho_c, rho_l=rho_l, rho_o=rho_o, tau_parm=tau_parm, kappa_c=kappa_c, kappa_l=kappa_l)
    aggregates = data.frame(I=I, c_S=c_S, c_I=c_I, c_R=c_R, l_S=l_S, l_I=l_I, l_R=l_R)
    tau_jt = tau_t("S", parms, aggregates)

    # print(aggregates)
    # print(tau_jt)

    value_S = total_utility(c_S,l_S,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha_U,Lbar) + discount_factor*((1-tau_jt)*continuation_value(SIR_, grid_list, contval_list$S) + tau_jt*continuation_value(SIR_, grid_list, contval_list$I))

    value_I = total_utility(c_I,l_I,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha_U,Lbar) + discount_factor*((1-pi_r-pi_d)*continuation_value(SIR_, grid_list, contval_list$I) + pi_r*continuation_value(SIR_, grid_list, contval_list$R)) + pi_d*U_D

    value_R = total_utility(c_R,l_R,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha_U,Lbar) + discount_factor*(continuation_value(SIR_, grid_list, contval_list$R))

    if(plannertype=="unweighted") {value = value_S + value_I + value_R}
    if(plannertype=="weighted") {value = S*value_S + I*value_I + R*value_R}

    return( value )
}

# interpolates choices for each type given: the agent type, the state, the state grid as a list, a dataframe of solved choices on the state grid
project_choices <- function(agent_type, state, grid_list, choices) {

    S = state$S
    I = state$I
    R = state$R

    if(agent_type=="S") {
        state_S = c(S,I,R)
        choices_S = choices$l_S

        l_S_fn = ipol(choices_S, grid=grid_list, method="multilinear")
        choice = l_S_fn(state_S)
    }
    if(agent_type=="I") {
        state_I = c(S,I,R)
        choices_I = choices$l_I

        l_I_fn = ipol(choices_I, grid=grid_list, method="multilinear")
        choice = l_I_fn(state_I)
    }
    if(agent_type=="R") {
        state_R = c(S,I,R)
        choices_R = choices$l_R

        l_R_fn = ipol(choices_R, grid=grid_list, method="multilinear")
        choice = l_R_fn(state_R)
    }

    return(choice)
}

# function to redo a grid to be Chebyshev nodes
make_cheby <- function(input) {
    a <- min(input)
    b <- max(input)
    n <- length(input)
    cheby_nodes <- rep(-1,length=n)
    for(k in 1:n) {
        x <- k/n - 1/(2*n)
        y <- pi/(2*n)
        cheby_nodes[k] <- 0.5*(a+b) + 0.5*(b-a)*sec(y)*cospi(x)
    }
    cheby_nodes <- sort(cheby_nodes)
    return(cheby_nodes)
}

# helper to build a uniform grid extension between two initial points. returns the (extension \setminus initial_points), which is of the chosen length
build_grid_extension <- function(startpoint, endpoint, length) {
    closed_length = length + 2
    extension = seq(from=startpoint, to=endpoint, length.out=closed_length)
    extension_interior =  setdiff(extension,c(startpoint,endpoint))
    return(extension_interior)
}

# builds a rectangular grid
build_grid <- function(Sgridlength, Igridlength, Rgridlength, cheby=0) {
    I_ext_threshold = 21
    I_ext_threshold_2 = 25

    S_ext_threshold = 21
    R_ext_threshold = 21

    S_base_piece <- seq(from=0, to=1, length.out=Sgridlength)
    I_base_piece <- seq(from=0, to=1, length.out=Igridlength)
    R_base_piece <- seq(from=0, to=1, length.out=Rgridlength)
    ifelse(cheby==1, S_base_piece <- make_cheby(S_base_piece), S_base_piece <- S_base_piece)
    ifelse(cheby==1, I_base_piece <- make_cheby(I_base_piece), I_base_piece <- I_base_piece)
    ifelse(cheby==1, R_base_piece <- make_cheby(R_base_piece), R_base_piece <- R_base_piece)

    grid <- as.data.frame(expand.grid(S_base_piece,I_base_piece,R_base_piece))
    
    colnames(grid) <- c("S","I","R")

    grid_list <- list(S = S_base_piece, I = I_base_piece, R = R_base_piece, dfrm = grid)
    return(grid_list)
}

# wrapper for SIR_series. generates a time series of SIR outcomes from a dataframe of policy functions
generate_time_series <- function(solved_values, interpolator="multilinear", grid_dfrm=NULL,...) {
    others_choices = data.frame(l_S = solved_values$labor_supply_S, l_I = solved_values$labor_supply_I, l_R = solved_values$labor_supply_R, simplex_check=simplex_check)

    others_utilities = data.frame(u_S = solved_values$lifetime_utility_S, l_I = solved_values$lifetime_utility_I, l_R = solved_values$lifetime_utility_R)

    choices_list = list(grid_list = grid_list, choices = others_choices, utilities = others_utilities)
    projection = SIR_series(exog_parms, choices_list, SIR_init, method=interpolator, grid_dfrm=grid_dfrm)

    results = cbind(time=c(1:final.time),projection)

    return(results)
}

# Generates a simple time series from an initial point, parameterization, and interpolator. Should output SIR, labor supply, consumption, and contact series.
simple_SIR_series <- function(parms,initial_condition,interpolators,sim_time) {
    w = parms$w
    phi = parms$phi
    final.time = parms$final.time
    pi_r = parms$pi_r
    pi_d = parms$pi_d
    rho_c = parms$rho_c
    rho_l = parms$rho_l
    rho_o = parms$rho_o

    S_0 = initial_condition[1]
    I_0 = initial_condition[2]
    R_0 = initial_condition[3]
    D_0 = 1 - (S_0 + I_0 + R_0)

    # state = c(S=S_0,I=I_0,R=R_0,benchmark_cons_contacts=impulse_init)
    state = data.frame(S=S_0,I=I_0,R=R_0,D=D_0) #for contact() and compute_R0()
    state_small = c(S=state[[1]],I=state[[2]],R=state[[3]]) #for interpolators
    types = c("S","I","R")
    labor_supply = data.frame(S=NA,I=NA,R=NA)
    # lifetime_utility = data.frame(S=NA,I=NA,R=NA)

    l_S = interpolators[["ls_fn"]]
    l_I = interpolators[["li_fn"]]
    l_R = interpolators[["lr_fn"]]

    l_S0 = l_S(state_small)
    l_I0 = l_I(state_small)
    l_R0 = l_R(state_small)

    choices_0 = data.frame(c_S = l_S(state_small)*w, c_I = l_I(state_small)*w*phi, l_S = l_S(state_small), l_I = l_I(state_small))

    initial_contacts = contacts(parms, choices_0)
    R0 = compute_R0(parms, choices_0, state)

    output_series = matrix(-1,nrow=sim_time, ncol=14)
    output_series[1,] = c(S_0,I_0,R_0,D_0,l_S0,l_I0,l_R0,w*l_S0,w*phi*l_I0,w*l_R0,(l_S0*S_0 + l_I0*I_0 + l_R0*R_0),(w*l_S0*S_0 + w*phi*l_I0*I_0 + w*l_R0*R_0),initial_contacts,R0)

    # output_series = matrix(rep(S_0,length.out=sim_time),
    #                             rep(I_0,length.out=sim_time),
    #                             rep(R_0,length.out=sim_time),
    #                             rep(D_0,length.out=sim_time),
    #                             rep(l_S0,length.out=sim_time),
    #                             rep(l_I0,length.out=sim_time),
    #                             rep(l_R0,length.out=sim_time),
    #                             rep(w*l_S0,length.out=sim_time),
    #                             rep(w*phi*l_I0,length.out=sim_time),
    #                             rep(w*l_R0,length.out=sim_time),
    #                             rep((l_S0*S_0 + l_I0*I_0 + l_R0*R_0),length.out=sim_time),
    #                             rep((w*l_S0*S_0 + w*phi*l_I0*I_0 + w*l_R0*R_0),length.out=sim_time),
    #                             rep(initial_contacts,length.out=sim_time),
    #                             rep(R0,length.out=sim_time), byrow=TRUE) 
    output_series <- as.data.frame(output_series)
    colnames(output_series) = c("S","I","R","D","labor_S","labor_I","labor_R","consumption_S", "consumption_I", "consumption_R", "aggregate_labor_supply", "aggregate_consumption", "contacts", "Reff")

    state = SIR_update(parms, choices_0, state)

    for(i in 2:sim_time) {
        state_small = c(S=state[[1]],I=state[[2]],R=state[[3]])
        labor_supply$S = l_S(state_small)
        labor_supply$I = l_I(state_small)
        labor_supply$R = l_R(state_small)

        choices = data.frame(c_S = labor_supply$S*w, c_I = w*phi*labor_supply$I, l_S = labor_supply$S, l_I = labor_supply$I)
        current_contacts = contacts(parms, choices)
        aggregate_labor_supply = labor_supply$S*state_small[[1]] + labor_supply$I*state_small[[2]] + labor_supply$R*state_small[[3]]
        aggregate_consumption = w*labor_supply$S*state_small[[1]] + w*phi*labor_supply$I*state_small[[2]] + w*labor_supply$R*state_small[[3]]

        new_row = c(S = state$S,
                                I = state_small[[1]],
                                R = state_small[[2]],
                                D = state_small[[3]],
                                labor_S = labor_supply$S,
                                labor_I = labor_supply$I,
                                labor_R = labor_supply$R,
                                consumption_S = labor_supply$S*w,
                                consumption_I = labor_supply$I*w*phi,
                                consumption_R = labor_supply$R*w,
                                aggregate_labor_supply = aggregate_labor_supply,
                                aggregate_consumption = aggregate_consumption,
                                contacts = current_contacts,
                                Reff = compute_R0(parms,choices,state))

        output_series[i,] = new_row

        state = SIR_update(parms, choices, state)
    }

    return(output_series)
}

# Generates impulse response functions
generate_IRF <- function(SIR_init, impulse_1_list, impulse_2_list, grid_list=NULL, impulse_time, final_time) {
    
    exog_parms_1 = impulse_1_list$exog_parms
    series_1 = simple_SIR_series(exog_parms_1, SIR_init, impulse_1_list, (impulse_time-1))

    print(head(series_1))

    new_init_cond = series_1[(impulse_time-1), 1:4]

    exog_parms_2 = impulse_2_list$exog_parms
    series_2 = simple_SIR_series(exog_parms_2, SIR_init, impulse_2_list, (final_time - impulse_time + 1))

    series = rbind(series_1,series_2)
    results = cbind(time=seq(1:final_time),series)

    return(results)
}

# wrapper for SIR_series. generates a time series of SIR outcomes from a dataframe of policy functions
generate_time_series_policy <- function(solved_values, policy_list, ...) {
    others_choices = data.frame(l_S = solved_values$labor_supply_S, l_I = solved_values$labor_supply_I, l_R = solved_values$labor_supply_R)

    others_utilities = data.frame(u_S = solved_values$lifetime_utility_S, l_I = solved_values$lifetime_utility_I, l_R = solved_values$lifetime_utility_R)

    choices_list = list(grid_list = grid_list, choices = others_choices, utilities = others_utilities)
    projection = SIR_series_policy(exog_parms, choices_list, SIR_init, policy_list)

    results = cbind(time=c(1:final.time),projection)

    return(results)
}


# runs the dynamic programming solver loop
dp_solver <- function(...) {

    ## initialize solver deltas and set tolerances
    damping_factor = 1

    VFI_count = 1

    l_lowerbound = 1e-6
    l_upperbound = Lbar*0.99
    l_upperbound_uniroot = Lbar*0.99
    epsilon_VFI = 0.005*mean(as.numeric(contval_list$S))
    if(problemtype=="planner"){epsilon_VFI = 0.0001*mean(as.numeric(contval_list$SWF))}
    #if(problemtype=="planner"){epsilon_VFI = 0.005*mean(as.numeric(contval_list$SWF))} #only for TESTING; tolerance is too low, leaves artifacts in S-type policy functions
    if(hot_start==1&problemtype=="planner") {epsilon_VFI = 0.05*epsilon_VFI}
    if(precision>0) {epsilon_VFI = (0.1^precision)*epsilon_VFI}

    delta_VFI = 100*epsilon_VFI
    delta_VFI_old = 2*delta_VFI
    simplex_check_1 = which(simplex_check==1)

    message("Labor supply choice bounded in (", l_lowerbound,", ", l_upperbound,"). VFI stopping threshold is ",epsilon_VFI,".")

    # begin outermost loop
            
    #sink("log.VFI.txt")
    while(delta_VFI > epsilon_VFI) {

        VFI_count_time = proc.time()[3]
        message("Beginning value function iteration ", VFI_count, " with ", sum(simplex_check), " simplex nodes")

        old_contval_mat = matrix( c(as.numeric(contval_list$S), as.numeric(contval_list$I), as.numeric(contval_list$R)), byrow=FALSE, ncol=3)

        old_planval_mat = matrix( c( as.numeric(contval_list$SWF) ), byrow=FALSE, ncol=1)
            
            #if(planner_solve==0) {
            VFI_results <- foreach(i=which(simplex_check==1), .export=ls(globalenv()), .packages=c("chebpol"), .combine=rbind, .inorder=TRUE) %dopar% {

                # begin innermost loop:
                delta_consistency = 10
                epsilon_consistency = 0.1

                state_entry_i = grid_dfrm[i,]
                others_choices_i = others_choices[i,]
                labor_supplies_new = c(others_choices_i[1], others_choices_i[2], others_choices_i[3])
                utilities_new = c(guess_dfrm$S_utility[i],guess_dfrm$I_utility[i],guess_dfrm$R_utility[i])

                iteration_consistency = 1
                if(problemtype=="planner") {
                        labor_supplies_old = others_choices_i
                        output = optim(par = as.numeric(labor_supplies_old), fn = pre_value_fn_planner, exog_parms = exog_parms, SIR_state = state_entry_i, grid_list = grid_list, contval_list = contval_list, control = list(fnscale=-1), method = "L-BFGS-B", lower=l_lowerbound, upper=l_upperbound)
            
                        print(output)

                        l_new = output$par

                        u_new = output$value
                        others_choices_i = l_new
                        print(l_new)

                        labor_supplies_new = l_new
                        print(labor_supplies_new)

                        others_choices_i = data.frame(S=labor_supplies_new[1], I=labor_supplies_new[2], R=labor_supplies_new[3])
            
                        value_S = pre_value_fn(choice = l_new[1], exog_parms = exog_parms, others_choices = others_choices_i, SIR_state = state_entry_i, agent_type = "S", grid_list = grid_list, contval_list = contval_list)
                        value_I = pre_value_fn(choice = l_new[2], exog_parms = exog_parms, others_choices = others_choices_i, SIR_state = state_entry_i, agent_type = "I", grid_list = grid_list, contval_list = contval_list)
                        value_R = pre_value_fn(choice = l_new[3], exog_parms = exog_parms, others_choices = others_choices_i, SIR_state = state_entry_i, agent_type = "R", grid_list = grid_list, contval_list = contval_list)

                        utilities_new = c(value_S, value_I, value_R, u_new)
            
                        return(c(labor_supplies_new, utilities_new))
                    }
                if(problemtype=="decentralized") {
                    while(delta_consistency > epsilon_consistency) {
    
                        if(iteration_consistency>1000) {damping_factor=0.5}
                        if(iteration_consistency>1500) {damping_factor=max(0.5^(iteration_consistency-1500),0.1) }
    
                        labor_supplies_old = others_choices_i

                        # calculate optimal labor supply for R-type
                        output_R = optim(par = labor_supplies_old[3], fn = pre_value_fn, exog_parms = exog_parms, others_choices = others_choices_i, SIR_state = state_entry_i, agent_type = "R", grid_list = grid_list, contval_list = contval_list, control = list(fnscale=-1), method = "L-BFGS-B", lower=l_lowerbound, upper=l_upperbound)
                        # update entries for R-type
                        l_R_new = as.numeric(damping_factor*output_R$par + (1-damping_factor)*labor_supplies_old[3])
                        u_R_new = as.numeric(damping_factor*output_R$value + (1-damping_factor)*utilities_new[3])
                        others_choices_i$R = l_R_new
    
                        # calculate optimal labor supply for I-type
                        output_I = optim(par = labor_supplies_old[2], fn = pre_value_fn, exog_parms = exog_parms, others_choices = others_choices_i, SIR_state = state_entry_i, agent_type = "I", grid_list = grid_list, contval_list = contval_list, control = list(fnscale=-1), method = "L-BFGS-B", lower=l_lowerbound, upper=l_upperbound)    
                        # update entries for S-type
                        l_I_new = as.numeric(damping_factor*output_I$par + (1-damping_factor)*labor_supplies_old[2])
                        u_I_new = as.numeric(damping_factor*output_I$value + (1-damping_factor)*utilities_new[2])
                        others_choices_i$I = l_I_new

                        if(measurezero==1) {
                            root_S = tryCatch(uniroot(m0_S, lower=0.01, upper=l_upperbound_uniroot, exog_parms = exog_parms, others_choices = others_choices_i, SIR_state = state_entry_i, grid_list = grid_list, contval_list = contval_list), error = function(err) {
                                object = data.frame(root=0)
                                return(object)
                                })
                            u_S_new = pre_value_fn(choice = root_S$root, exog_parms = exog_parms, others_choices = others_choices_i, SIR_state = state_entry_i, agent_type = "S", grid_list = grid_list, contval_list = contval_list)
                            l_S_new = as.numeric(damping_factor*root_S$root + (1-damping_factor)*labor_supplies_old[1])
                            u_S_new = as.numeric(damping_factor*u_S_new + (1-damping_factor)*utilities_new[1])
                        }
                        
                        # Update the labor supplies with the new choices and lifetime utilities with the new values
                        labor_supplies_new = c(l_S_new, l_I_new, l_R_new)
    
                        SWF_new = u_S_new + u_I_new + u_R_new

                        utilities_new = c(u_S_new, u_I_new, u_R_new, SWF_new)
                        others_choices_i = data.frame(S=labor_supplies_new[1], I=labor_supplies_new[2], R=labor_supplies_new[3])
    
                        delta_consistency = max(abs(labor_supplies_new - labor_supplies_old))
                        message("Finished consistency iteration ", iteration_consistency, ". Delta is ", round(delta_consistency,5), ", new labor supplies are (l_S, l_I, l_R) = (", round(l_S_new,3), ", ", round(l_I_new,3), ", ", round(l_R_new,3), "). \n")
                        iteration_consistency = iteration_consistency + 1
                #sink()
                    }
                return(c(labor_supplies_new, utilities_new))
                }
            }

        VFI_count_time = proc.time()[3] - VFI_count_time
        if(problemtype=="decentralized"){
            new_contval_mat = old_contval_mat
    
            new_contval_mat[simplex_check_1,] = c(as.numeric(VFI_results[,4]), as.numeric(VFI_results[,5]), as.numeric(VFI_results[,6]) )
    
            delta_VFI_old = delta_VFI
            delta_VFI = max(abs(new_contval_mat - old_contval_mat))
        }
        if(problemtype=="planner") {
            new_planval_mat = old_planval_mat

            new_planval_mat[simplex_check_1] = c(as.numeric(VFI_results[,7]))

            delta_VFI_old = delta_VFI
            delta_VFI = max(abs(new_planval_mat - old_planval_mat))
            if( (delta_VFI_old - delta_VFI < min(1e-04,epsilon_VFI) ) & (delta_VFI_old - delta_VFI >= 0) & (VFI_count > 100) ) { delta_VFI = epsilon_VFI/2}
            #if( (VFI_count > 2) ) { delta_VFI = 1e-04}
        }
        message("Finished iteration ", VFI_count, ". Delta is ", round(delta_VFI,10) ,". Change in delta is ", (delta_VFI_old - delta_VFI) ,". Total time: ", round(VFI_count_time,3) ," seconds. Average core-time per state: ", ncores*round(VFI_count_time/nrow(grid_dfrm),3), " seconds. Average value function levels: (S,I,R) = (", round(mean(VFI_results[,4]),1),",", round(mean(VFI_results[,5]),1), ",", round(mean(VFI_results[,6]),1),"), SWF = ", round(mean(VFI_results[,7]),1),".") 

        intermediate_results <- as.data.frame(VFI_results)
        colnames(intermediate_results) <- c("labor_supply_S","labor_supply_I","labor_supply_R","lifetime_utility_S","lifetime_utility_I","lifetime_utility_R","SWF")

        if(precision>1) {fwrite(intermediate_results, file=paste0("../../Results/value_policy_functions/value_policy_functions__",scenario_label,".csv"))}

        assign.time = proc.time()[3]
        message("Assigning vectors...") 
        contval_list$S[simplex_check_1] = as.numeric(VFI_results[,4])
        contval_list$I[simplex_check_1] = as.numeric(VFI_results[,5])
        contval_list$R[simplex_check_1] = as.numeric(VFI_results[,6])
        contval_list$SWF[simplex_check_1] = as.numeric(VFI_results[,7])

        VFI_count = VFI_count + 1

        assign.time = round(proc.time()[3] - assign.time,3)

        message("Finished assigning vectors. Time taken: ", assign.time, " seconds.")

    }

    # create "solved_values", initialized as the starting values
    solved_values = data.frame(labor_supply_S=others_choices_init$S, labor_supply_I=others_choices_init$I, labor_supply_R=others_choices_init$R, lifetime_utility_S=contval_list_init$S, lifetime_utility_I=contval_list_init$I, lifetime_utility_R=contval_list_init$R, SWF = contval_list$SWF, simplex_check=simplex_check)

# print(dim(VFI_results))
    # update "solved_values" with VFI_results
    VFI_results_numeric = cbind(matrix(as.numeric(VFI_results),byrow=FALSE,ncol=7),simplex_check[simplex_check_1])
# print(dim(VFI_results_numeric))

    solved_values[simplex_check_1,] = VFI_results_numeric

    return(solved_values)

}

# wrapper to generate a list of plots from a time series of outcomes
generate_plots <- function(results) { 
    time_series_plots_base = ggplot(results, aes(x=time))

    aggregate_SRs = time_series_plots_base + 
                    geom_line(aes(y=S), color="black", size=1) +
                    geom_line(aes(y=R), color="dodgerblue3", size=1) +
                    theme_bw() + ggtitle("SR sizes (S,R) = (black,blue)") +
                    xlab("Day") + ylab("Proportion of population") +
                    ylim(c(0,1))

    aggregate_IDs = time_series_plots_base + 
                    geom_line(aes(y=I), color="firebrick4", size=1) +
                    geom_line(aes(y=D), color="darkgray", size=1) +
                    theme_bw() + ggtitle("ID sizes (I,D) = (red,gray)") +
                    xlab("Day") + ylab("Proportion of population")

    aggregate_deviations = time_series_plots_base + 
                    geom_line(aes(y=aggregate_hours_deviation), color="black", size=1) +
                    geom_line(aes(y=aggregate_consumption_deviation), color="dodgerblue4", size=1) +
                    theme_bw() + ggtitle("Aggregate hours (black) and consumption (blue)") +
                    xlab("Day") + ylab("% deviations from initial steady state")

    individual_l = time_series_plots_base + 
                    geom_line(aes(y=labor_S), color="black", size=1) +
                    geom_line(aes(y=labor_I), color="firebrick4", size=1) +
                    geom_line(aes(y=labor_R), color="dodgerblue3", size=1) +
                    theme_bw() + ggtitle(paste0("Labor supplies. Contact utility weights (c,l) = ", round(gamma_c,4), ", ", round(gamma_l,4),", CRRA = ", risk_aversion,".") ) +
                    xlab("Day") + ylab("Proportion of available time") # +
                   # ylim(c(0,0.25))

    Reff_plot = time_series_plots_base +
                    geom_line(aes(y=Reff), color="black", size=1) +
                    geom_hline(yintercept = 2.6, linetype="dashed", color="darkgray") +
                    theme_bw() + ggtitle("R_effective") + 
                    xlab("Day") + ylab("Number infected")

    total_cases_plot = time_series_plots_base +
                    geom_line(aes(y=total_cases), color="black", size=1) +
                    theme_bw() + ggtitle("Cumulative cases till date") +
                    xlab("Day") + ylab("Proportion of population")

    infection_prob_plot = time_series_plots_base +
                    geom_line(aes(y=prob_infection_S), color="black", size=1) +
                    theme_bw() + ggtitle("Probability of getting infected") + 
                    xlab("Day") + ylab("Probability")                    

    weekly_new_cases_plot = time_series_plots_base +
                    geom_line(aes(y=weekly_new_cases), color="black", size=1) +
                    theme_bw() + ggtitle("Daily new cases") +
                    xlab("Day") + ylab("Proportion of population")

    pop_weighted_labor_supplies = time_series_plots_base +
                    geom_line(aes(y=labor_S*S), color="black", size=1) +
                    geom_line(aes(y=labor_I*I), color="firebrick4", size=1) +
                    geom_line(aes(y=labor_R*R), color="dodgerblue3", size=1) +
                    theme_bw() + ggtitle("Population-weighted labor supplies") +
                    xlab("Day") + ylab("Hours worked by type")

    total_contacts = time_series_plots_base + 
                    geom_line(aes(y=(consumption_contacts+labor_contacts+other_contacts)), color="black", size=1) +
                    geom_hline(yintercept = 12, linetype="dashed", color="darkgray") +
                    theme_bw() + ggtitle("Total contacts") +
                    xlab("Day") + ylab("Daily contacts")

    contacts_by_activity = time_series_plots_base +
                    geom_line(aes(y=consumption_contacts), color="#1b9e77", size=1) +
                    geom_line(aes(y=labor_contacts), color="#d95f02", size=1) +
                    geom_line(aes(y=other_contacts), color="#7570b3", size=1) +
                    theme_bw() + ggtitle("Contacts at consumption, labor, other (green, orange, purple)") +
                    xlab("Day") + ylab("Contacts by activity")

    #list_of_plots = list(aggregate_SIRs, Reff_plot, individual_l, weekly_new_cases_plot, aggregate_deviations, total_cases_plot, pop_weighted_labor_supplies, contacts_by_activity)
  
    list_of_plots = list(aggregate_SRs, aggregate_IDs, Reff_plot, 
                        individual_l, aggregate_deviations, infection_prob_plot, 
                        pop_weighted_labor_supplies, contacts_by_activity, weekly_new_cases_plot)

    return(list_of_plots)
}

# expands (or shrinks) a solved value function using linear interpolation on the smaller grid to project to the larger grid
grid_grow <- function(solved_values, new_grid_dfrm, simplex_check) {
    value_fn_S_old = solved_values$lifetime_utility_S
    value_fn_I_old = solved_values$lifetime_utility_I
    value_fn_R_old = solved_values$lifetime_utility_R
    policy_fn_S_old = solved_values$labor_supply_S
    policy_fn_I_old = solved_values$labor_supply_I
    policy_fn_R_old = solved_values$labor_supply_R

    old_grid_list = list(S=unique(solved_values$S), I=unique(solved_values$I), R=unique(solved_values$R))

    new_data = data.frame(new_grid_dfrm, simplex_check = simplex_check)

    vfn_S_o = ipol(value_fn_S_old, grid=old_grid_list, method="multilinear")
    vfn_I_o = ipol(value_fn_I_old, grid=old_grid_list, method="multilinear")
    vfn_R_o = ipol(value_fn_R_old, grid=old_grid_list, method="multilinear")
    pfn_S_o = ipol(policy_fn_S_old, grid=old_grid_list, method="multilinear")
    pfn_I_o = ipol(policy_fn_I_old, grid=old_grid_list, method="multilinear")
    pfn_R_o = ipol(policy_fn_R_old, grid=old_grid_list, method="multilinear")

    colnames(new_grid_dfrm) = NULL
    new_grid_dfrm = t(as.matrix(new_grid_dfrm))

    vfn_S_n = vfn_S_o(new_grid_dfrm)
    vfn_I_n = vfn_I_o(new_grid_dfrm)
    vfn_R_n = vfn_R_o(new_grid_dfrm)
    pfn_S_n = pfn_S_o(new_grid_dfrm)
    pfn_I_n = pfn_I_o(new_grid_dfrm)
    pfn_R_n = pfn_R_o(new_grid_dfrm)

    new_values = data.frame(lifetime_utility_S = vfn_S_n, lifetime_utility_I = vfn_I_n, lifetime_utility_R = vfn_R_n, labor_supply_S=pfn_S_n, labor_supply_I=pfn_I_n, labor_supply_R=pfn_R_n, new_data)

    return(new_values)
}


##################################################################################################
################################# HELPER GRAPHS & FUNCTIONS ###################################### ----
##################################################################################################


# helper function to plot solved value and policy functions
plot_pfn_vfn <- function(vfn,pfn,basegrid,labels) {
    fv_mat <- t(matrix(vfn,nrow=length(basegrid)))
    po_mat <- t(matrix( pfn,nrow=length(basegrid)))

    image2D(z=fv_mat,x=basegrid,y=basegrid,xlab=c("Susceptibles"),ylab=c("Infecteds"),col=plasma(n=100),contour=TRUE,main=labels[1])
    image2D(z=po_mat,x=basegrid,y=basegrid,xlab=c("Susceptibles"),ylab=c("Infecteds"),col=plasma(n=100), main=labels[2])
}



# creates charts to assess value and policy functions
chart_value_fn <- function(df=solved_values, quantile_no=20,use_quantiles=1,no_quantile_start=0, no_quantile_end=1,type="S",var="lifetime_utility", title="plot",colors=1) {
    #creates dataframe of quantile information (quantiles) so that these can be used as thresholds to be plotted as separate series
    count=1
    quantiles = data.frame(count=0,Decile=0,Val=0,Other=0)
    increment = 1/quantile_no
    quant_range = seq(0,no_quantile_end,increment)
    for (val in quant_range)
    {
    dec = round(quantile(df[,paste0(var,"_", type)], probs = c(val)),digits=10)
    other = no_quantile_start + val
    quantiles[count,] = c(count,val,dec,other)
    count = count + 1
    }
    message(quantiles)
    
    #separates the state space points into different quantiles of the variable so these quantile ranges can be plotted as different colors, values of S are stored in x 
    #and values of I are stored in y, values of R are imputed as 1-S-I = 1-x-y
    #Note if use_quantiles set to 0 then just uses even spacing on (0,1) rather than interquartile ranges.
    x=list()
    y=list()
    quant_range = c(1:quantile_no)
    for(val in quant_range)
    {
    z <- quantiles$Val[which(quantiles$count==val)]
    if(use_quantiles==0) {  z <- quantiles$Other[which(quantiles$count==val)]    }
    x[val] <- list(df[,paste0("S")][which((df[,paste0(var,"_", type)]>z)&(df[,paste0("simplex_check")]==1))]) #this is a list of lists, one list for each quantile
    y[val] <- list(df[,paste0("I")][which((df[,paste0(var,"_", type)]>z)&(df[,paste0("simplex_check")]==1))]) #as prev line
    }


    #runs the simpex plot command using the two lists of lists, number of quantiles and labels
    value_plot = simplex(x,y,quantile_no,quantiles,var,use_quantiles,colors,label= expression(italic("S=1"), italic("I=1"), italic("R=1")), title=title)

    return(value_plot)
}


# function to plot value and policy functions on the S I R simplex
simplex <- function(x, y, quantile_no,quantiles, var, use_quantiles, colors,
                    label = expression(italic(x)[1], italic(x)[2], italic(x)[3]), title="title") {

    # creates empty square plot from two cartesian points (-0.2,-0.2) all the way to (1.2,1.2)
    #op <- par(mar=c(1,0,0,0) + 0.1, pty="s")
    op <- par(mar=c(1,1,1,1) + 0, pty="s")
    plot(x=c(-0.05,1.05), y=c(-0.05,1.05), type="n", axes=FALSE, xlab="", ylab="", main=title)

    # transform the points (i.e. a value for S and a value for I -- the value for R is then imputed as 1-S-I) stored in each list within the list of lists 
    # into simplex terms
    xx=x
    yy=y
    quant_range = c(1:quantile_no)
    for(val in quant_range)
    {
        xx[[val]] <- 0.5*(1-x[[val]]+y[[val]])
        yy[[val]] <- 0.5*sqrt(3)*(1-x[[val]]-y[[val]])
    }

    # generate list of colors -- two groups, first half of quantiles (0% to 50%) are in color_list1 which goes from dark blue (cold) to white (middle); second half (50%
    # -100%) are in color_list2 which goes from white to orange (hot). I.e. 50th percentile represented by white points.
    colfunc <- colorRampPalette(c("#3936C9", "#FFFFFF"))  #"#E4E4F5"
    color_list1 = colfunc(quantile_no/2)
    colfunc <- colorRampPalette(c("#FFFFFF", "#E88813"))
    color_list2 = colfunc(quantile_no/2)
    color_list = c(color_list1,"#FFFFFF",color_list2)
    
    #or simple color list if colors set to 0
    if(colors==0) {
      colfunc <- colorRampPalette(c("#FFFFFF", "#000000"))
      color_list = colfunc(quantile_no)
    }

    # plot points -- points are selected to be large enough (cex=4) so that they cover all remaining spaces creating the look of an area graph, higher quantiles are
    # layered over the top of lower quantiles.
    for(val in quant_range) 
    {
        points(x=xx[[val]], y=yy[[val]], type='p',pch=16,cex=4,col=color_list[val])
    }

    #get rid of excess -- due to size of points, the points on the border of the simplex will go over outside of the simplex. Use these polygons to cover this space white.
    #polygon axes are chosen to create 4 sided shapes running along the axes -- note x represents x coordinate of the four points making up the four-sided polygon, y
    #represents the y coordinates of these four points. The lines created in the next set of the code describe the triangle borders and hence were used to find the 
    #coordinates of the polygons here.
    polygon(x = c(0,0.5,-0.2,-0.2), y = c(0,0.5*sqrt(3),1.2,0),col ="#FFFFFF",border="#FFFFFF")
    polygon(x = c(1,0.5,1.2,1.2), y = c(0,0.5*sqrt(3),1.2,0),col ="#FFFFFF",border="#FFFFFF")
    polygon(x = c(-0.2,1.2,-0.2,1.2), y = c(0,0,-0.5,-0.5),col ="#FFFFFF",border="#FFFFFF")
    polygon(x = c(0.3,0.7,0.3,0.7), y = c(0.5*sqrt(3),0.5*sqrt(3),0.5*sqrt(3)+0.5,0.5*sqrt(3)+0.5),col ="#FFFFFF",border="#FFFFFF")

    # triangle (borders) -- cartesian points of triangle that draw out simplex are (0,0) -leftmost point; (1,0)-rightmost point; (0.5,0.5*sqrt(3))-top point.
    points(x=c(0,0.5,1,0), y=c(0,0.5*sqrt(3),0,0), type="l",lwd=2)

    #axis lines and labels
    ux = seq(0,1,0.2) #S = 0 to 1 in 0.2 increments
    uy = c(rep(0,6))  #I is 0
    vx = seq(0,1,0.2) #S = 0 to 1 in 0.2 increments
    vy = rev(vx)      #I=1-S such that Z=0

    uxx=ux
    uyy=uy
    vxx=vx
    vyy=vy

    #first converts points from SIR space into cartesian coordinates, then draws lines between each set of two points, then adds labels to each axis
    for(val in c(1:6))
    {
    uxx[val] <- 0.5*(1-ux[val]+uy[val])
    uyy[val] <- 0.5*sqrt(3)*(1-ux[val]-uy[val])
    vxx[val] <- 0.5*(1-vx[val]+vy[val])
    vyy[val] <- 0.5*sqrt(3)*(1-vx[val]-vy[val])
    points(x=c(uxx[val],vxx[val]), y=c(uyy[val],vyy[val]), type="l")
    if(val != 6) {legend(x=vxx[val]-0.05,y=vyy[val]-0.02,legend="",bty="n",title=paste0("S=",ux[val]),cex=0.5,pt.cex=0.7)}
    if(val != 6) {legend(x=uxx[val]-0.1,y=uyy[val]+0.04,legend="",bty="n",title=paste0("S=",ux[val]),cex=0.5,pt.cex=0.7)}
    }

    # labels
    if (!is.null(label)) {
        text(x=0.5, y=0.5*sqrt(3), pos=3, labels=label[3])
        text(x=0.0, y=0.0, pos=2, labels=label[1])
        text(x=1.0, y=0.0, pos=4, labels=label[2])
    }

    
    ### generate legend
    #for rounding
    i=2
    if(var=="lifetime_utility") {i=0}
  
    bottom_val = round(quantiles$Val[which(quantiles$count==1)],digits=i)
    middle_val = round(quantiles$Val[which(quantiles$count==11)],digits=i)
    top_val = round(quantiles$Val[which(quantiles$count==21)],digits=i)

    # if not using quantiles then straight percentile labels (stored in quantiles$Decile) are used e.g. 0.05 0.1 0.15 etc
    if(use_quantiles==0) {
      bottom_val = round(quantiles$Other[which(quantiles$count==1)],digits=i)
      middle_val = round(quantiles$Other[which(quantiles$count==11)],digits=i)
      top_val = round(quantiles$Other[which(quantiles$count==21)],digits=i)
    }
    
    
    legend_labels = c(bottom_val,rep(NA,8),middle_val,rep(NA,9),top_val)
    legend(x=0.8,y=0.5*sqrt(3),legend=rev(legend_labels),fill=rev(color_list),border="#FFFFFF",bty="n",y.intersp
    =0.3,cex=0.8,pt.cex=1)
    legend(x=0.8,y=0.5*sqrt(3)+0.08,legend="",bty="n",title=var,title.adj=1,cex=0.9,pt.cex=1)

    # restore plotting parameters
    par(op)
}

fig2_sketch <- function(plot_base, groupvar, summary_table) {

    infection = plot_base + 
                geom_line(aes(y=I, group=get(groupvar), color=get(groupvar)), size=1) +
                geom_line(aes(y=D, group=get(groupvar), color=get(groupvar)), linetype="dashed", size=1) +
                theme_bw() + ggtitle("Infecteds (solid) and dead (dashed)") +
                xlab("Day") + ylab("Proportion of population") + 
                scale_colour_manual(values=Mix) + guides(color=FALSE)

    Reff = plot_base + 
                geom_line(aes(y=Reff, group=get(groupvar), color=get(groupvar)), size=1) +
                theme_bw() + ggtitle("R_effective") +
                xlab("Day") + ylab("R_effective") + 
                geom_hline(yintercept=1, linetype="dashed", alpha=0.75, color="darkgray") +
                scale_colour_manual(values=Mix) + guides(color=FALSE)

    cases = plot_base + 
                geom_line(aes(y=weekly_new_cases, group=get(groupvar), color=get(groupvar)), size=1) +
                theme_bw() + ggtitle("Daily new cases") +
                xlab("Day") + ylab("Cases") + 
                scale_colour_manual(values=Mix) + guides(color=FALSE)

    infect_prob_S = plot_base + 
                geom_line(aes(y=prob_infection_S, group=get(groupvar), color=get(groupvar)), size=1) +
                theme_bw() + ggtitle("Probability of infection") +
                xlab("Day") + ylab("Probability") + 
                scale_colour_manual(values=Mix) + guides(color=FALSE)

    uninfected = plot_base + 
                geom_line(aes(y=(S+R), group=get(groupvar), color=get(groupvar)), size=1) +
                theme_bw() + ggtitle("Uninfecteds (effective labor force)") +
                xlab("Day") + ylab("Proportion of population") + 
                scale_colour_manual(values=Mix) + guides(color=FALSE)

    recovereds = plot_base + 
                geom_hline(yintercept = 1 - 1/all_orderings$Reff[1], linetype = "dashed") +
                geom_line(aes(y=R, group=get(groupvar), color=get(groupvar)), size=1) +
                theme_bw() + ggtitle("Total recovered") +
                xlab("Day") + ylab("Proportion of population") + 
                scale_colour_manual(values=Mix) + guides(color=FALSE)

    consumption = plot_base + 
                geom_line(aes(y=aggregate_consumption_deviation, group=get(groupvar), color=get(groupvar)), size=1) +
                theme_bw() + ggtitle("GDP deviation") +
                xlab("Day") + ylab("Percentage of initial steady state") + 
                scale_colour_manual(values=Mix) + labs(linetype="Contact\nordering")

    hours = plot_base + 
                geom_line(aes(y=aggregate_hours_deviation, group=get(groupvar), color=get(groupvar)), size=1) +
                theme_bw() + ggtitle("Aggregate hours deviation") +
                xlab("Day") + ylab("Percentage of initial steady state") + 
                scale_colour_manual(values=Mix)

    weighted_labor_supplies = plot_base +
                geom_line(aes(y = S*labor_S, group=get(groupvar), linetype=get(groupvar)), size=1, color = "black") +
                geom_line(aes(y = I*labor_I, group=get(groupvar), linetype=get(groupvar)), size=1, color = "firebrick4") +
                geom_line(aes(y = R*labor_R, group=get(groupvar), linetype=get(groupvar)), size=1, color = "dodgerblue4") +
                theme_bw() + ggtitle("Aggregate labor supplies\n(black = S, red = I, blue = R)") +
                xlab("Day") + ylab("pop_size*hours") + theme(legend.key.size=grid::unit(2,"lines")) + labs(linetype="Contact\nordering")

    total_contacts = plot_base + 
                geom_line(aes(y=total_contacts, group=get(groupvar), color=get(groupvar)), size=1) +
                geom_hline(yintercept = 12, linetype="dashed", color="darkgray") +
                geom_hline(yintercept = 5, linetype="dashed", color="darkgray") +
                theme_bw() + ggtitle("Total S-I contacts") +
                xlab("Day") + ylab("Daily contacts") + 
                scale_colour_manual(values=Mix) + guides(color=FALSE)

    labor_S = plot_base + 
                geom_line(aes(y=labor_S, group=get(groupvar), color=get(groupvar)), size=1) +
                theme_bw() + ggtitle("labor supply S") +
                xlab("Day") + ylab("Time spent") + 
                scale_colour_manual(values=Mix) + guides(color=FALSE)

    labor_I = plot_base + 
                geom_line(aes(y=labor_I, group=get(groupvar), color=get(groupvar)), size=1) +
                theme_bw() + ggtitle("labor supply I") +
                xlab("Day") + ylab("Time spent") + ylim(c(2,8)) +
                scale_colour_manual(values=Mix) + guides(color=FALSE)

    labor_R = plot_base + 
                geom_line(aes(y=labor_R, group=get(groupvar), color=get(groupvar)), size=1) +
                theme_bw() + ggtitle("labor supply R") +
                xlab("Day") + ylab("Time spent") + ylim(0,max(all_orderings$labor_R)) +
                scale_colour_manual(values=Mix) + guides(color=FALSE)

    list_of_plots = list(infection, recovereds, Reff, total_contacts, labor_S, labor_I, consumption, weighted_labor_supplies, summary_table)

    return(list_of_plots)

}
