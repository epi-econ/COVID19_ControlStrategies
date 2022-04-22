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
SIR_series <- function(parms, choices_list, aggregates, method="multilinear", grid_dfrm=NULL, grid_list=NULL, infolag_list=NULL, fallback_choices_list=NULL, policy_list=NULL, solved_values=NULL, fallback_solved_values=NULL) {
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
    test_quality <- 1
    if("test_quality"%in%names(parms)==TRUE) {test_quality <- parms$test_quality}
    tq_seq <- NULL
    policy_switch <- 0
    problem_type <- as.character(unique(solved_values$type))

    if("infolag"%in%names(parms)==TRUE) { infolag <- parms$infolag }
    if("infolag"%in%names(parms)==FALSE) { infolag <- 0 }
    if(length(infolag_list)>0) {
        infolag_list <- infolag_list
        infolag_vec <- infolag_list$lag_sequence
    }

    if("quality_sequence"%in%names(infolag_list)==TRUE) { tq_seq <- infolag_list$quality_sequence }

    if(is.null(fallback_choices_list)==FALSE) { 
        # message("yeet1")
        fallback_labor_supply <- data.frame(S=NA,I=NA,R=NA) 
        fallback_grid_list = fallback_choices_list$grid_list
        fallback_choices = fallback_choices_list$choices
    }

    if("compliance"%in%names(parms)==TRUE) { compliance <- parms$compliance }
    if("compliance"%in%names(parms)==FALSE) { compliance <- 1 }

    if("eqm_concept"%in%names(parms)==TRUE) { 
        ncores <- min(detectCores()-1,18)
        cl = makeCluster(ncores, type="FORK")
        registerDoParallel(cl)
        setDefaultCluster(cl=cl)

        message("yeetleets")

        eqm_concept <- parms$eqm_concept 
        mixing_prob_mat <- as.data.frame(matrix(NA,nrow=final.time,ncol=9))
        mixing_prob_mat[1,] <- c(1,0,0, 0,1,0, 0,0,1)
        colnames(mixing_prob_mat) <- c("l_SS","l_SI","l_SR","l_IS","l_II","l_IR","l_RS","l_RI","l_RR")
    }
    if("eqm_concept"%in%names(parms)==FALSE) { 
        eqm_concept <- "default" 
        mixing_prob_mat <- as.data.frame(matrix(NA,nrow=final.time,ncol=9))
        mixing_prob_mat[1,] <- c(1,0,0, 0,1,0, 0,0,1)
        colnames(mixing_prob_mat) <- c("l_SS","l_SI","l_SR","l_IS","l_II","l_IR","l_RS","l_RI","l_RR")
    }

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

    if(length(policy_list)>0) {
        policy_switch = policy_list$switch
        policy_start = policy_list$start_date
        policy_stop = policy_list$end_date
        policy_labor_restriction = policy_list$labor_supply_level
        policy_state_dependent = policy_list$state_dependent
        restricted_level = l_S0*policy_labor_restriction
    }

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

    perceived_state <- state
    mixing_probs <- c(1,0,0, 0,1,0, 0,0,1)

    qre.tm.vec <- rep(NA,length.out=final.time)

    for(i in 2:final.time) {
    perceived_state <- state

    if(length(infolag_list)>0) {infolag <- infolag_vec[i]}
        # if(method=="multilinear") {
            for(j in 1:length(types)) {
                choices = choices_list$choices
                agent_type = types[j]

                # apply information lag
                if(length(infolag_list)>0) {
                    if(i<=infolag&infolag>0) {
                        perceived_state <- data.frame(S=S_0,I=I_0,R=R_0)
                    }
                    if(i>infolag&infolag>0) {
                        perceived_state <- data.frame(S=output_series$S[(i-infolag)], I=output_series$I[(i-infolag)], R=output_series$R[(i-infolag)])
                        # print("dink")
                    }
                }
                if(("infolag"%in%names(parms))&(length(infolag_list)==0)) {
                    if(i<=infolag&infolag>0) {
                        perceived_state <- data.frame(S=S_0,I=I_0,R=R_0)
                    }
                    if(i>infolag&infolag>0) {
                        perceived_state <- data.frame(S=output_series$S[(i-infolag)], I=output_series$I[(i-infolag)], R=output_series$R[(i-infolag)])
                    }
                }

                # project choices
                labor_supply[j] = project_choices(agent_type, perceived_state, choices_list$grid_list, choices)

                if(is.null(fallback_choices_list)==FALSE) {
                    # message("yeet")
                    fallback_labor_supply[j] <- project_choices(agent_type, perceived_state, fallback_grid_list, fallback_choices)
                }

                if((length(policy_list>0))){
                    if((policy_switch==1)) {
                        if((i>=policy_start)&(i<=policy_stop)){
                            labor_supply[j] = min(labor_supply[j],restricted_level)
                        }
                    }
                }
            }
        # }
        
        # message("Period ", i,". True state is: ", paste(state,","))
        # message("Period ", i,". Perceived state is: ", paste(perceived_state,","))

        l_S = labor_supply$S
        l_I = labor_supply$I
        l_R <- labor_supply$R

        # Apply the test quality effect
        if(length(tq_seq>0)) {
            test_quality <- tq_seq[i]
        }
        # message("Period ", i,". Correct actions are: ", paste(labor_supply,","))
        if(eqm_concept!="qre"){
            l_S <- test_quality*l_S + (1-test_quality)*( c(1/3,1/3,1/3)%*%c(labor_supply$S, labor_supply$I, labor_supply$R) )
            l_I <- test_quality*l_I + (1-test_quality)*( c(1/3,1/3,1/3)%*%c(labor_supply$S, labor_supply$I, labor_supply$R) )
            l_R <- test_quality*l_R + (1-test_quality)*( c(1/3,1/3,1/3)%*%c(labor_supply$S, labor_supply$I, labor_supply$R) )
        }
        if(eqm_concept=="qre") {
            l_S <- labor_supply$S
            l_I <- labor_supply$I
            l_R <- labor_supply$R

            probs <- mixing_probs # start with everyone playing last period's mix
            choice_vector <- data.frame(S = l_S, I = l_I, R = l_R) # a choice vector of pure strats
            lambda <- 1/(1-test_quality) - 1

            message("Solving for QRE in period ",i,", lambda is ",round(lambda,2))
            # message(probs)
            # message(head(solved_values))
            # message(choice_vector)
            # message(perceived_state)
            qre.tm <- as.numeric(proc.time()[3])
            if(is.finite(lambda)){
                mixing_probs <- optimParallel(par=probs, fn=qre_system, method="L-BFGS-B", lower=0, upper=1, lambda=lambda, payoffs_SIR=solved_values, problem_type=problem_type, choices=choice_vector, SIR_state=perceived_state, exog_parms=exog_parms)$par
            }
            if(!is.finite(lambda)){
                mixing_probs <- c(1,0,0, 0,1,0, 0,0,1)
            } 
            qre.tm.vec[i] <- as.numeric(round(proc.time()[3] - qre.tm,2))
            message("QRE solved! Time taken: ", qre.tm.vec[i])

            mixing_prob_mat[i,] <- mixing_probs

            # print(as.matrix(mixing_probs[1:3],nrow=1))
            # print(as.matrix(choice_vector,ncol=1))
            # as.numeric(choice_vector)

            l_S <- sum(mixing_probs[1:3]*choice_vector)
            l_I <- sum(mixing_probs[4:6]*choice_vector)
            l_R <- sum(mixing_probs[7:9]*choice_vector)

            # print(l_S)
            # message("yeetfeet")
        }

        # Push the changes back to labor_supply
        labor_supply$S <- as.numeric(l_S)
        labor_supply$I <- as.numeric(l_I)
        labor_supply$R <- as.numeric(l_R)
        # message("Period ", i,". Chosen actions are: ", paste(labor_supply,","))

        # Apply the compliance effect
        if(is.null(fallback_choices_list)==FALSE) {
            # message("Period ", i,". Compliant actions are: ", paste(labor_supply,","))
            # message("yeet")

            fallback_l_S <- test_quality*fallback_labor_supply$S + (1-test_quality)*( c(1/3,1/3,1/3)%*%c(fallback_labor_supply$S, fallback_labor_supply$I, fallback_labor_supply$R) )
            fallback_l_I <- test_quality*fallback_labor_supply$I + (1-test_quality)*( c(1/3,1/3,1/3)%*%c(fallback_labor_supply$S, fallback_labor_supply$I, fallback_labor_supply$R) )
            fallback_l_R <- test_quality*fallback_labor_supply$R + (1-test_quality)*( c(1/3,1/3,1/3)%*%c(fallback_labor_supply$S, fallback_labor_supply$I, fallback_labor_supply$R) )
            # message("Non-compliant actions are: ", paste(c(fallback_l_S, fallback_l_I, fallback_l_R),","))

            labor_supply$S <- compliance*l_S + (1-compliance)*fallback_labor_supply$S
            labor_supply$I <- compliance*l_I + (1-compliance)*fallback_labor_supply$I
            labor_supply$R <- compliance*l_R + (1-compliance)*fallback_labor_supply$R
            # message("Period ", i,". Chosen actions are: ", paste(labor_supply,","))
            
            l_S <- labor_supply$S
            l_I <- labor_supply$I
            l_R <- labor_supply$R

            if(eqm_concept=="qre") {
                l_S <- fallback_labor_supply$S
                l_I <- fallback_labor_supply$I
                l_R <- fallback_labor_supply$R

                probs <- mixing_probs # start with everyone playing last period's mix
                choice_vector <- data.frame(S = l_S, I = l_I, R = l_R) # a choice vector of pure strats
                # lambda <- pmin(1/(1-test_quality) - 1, 1e10)
                lambda <- 1/(1-test_quality) - 1

                message("Solving for fallback QRE in period ",i,", lambda is ",lambda)
                qre.tm <- as.numeric(proc.time()[3])
                if(is.finite(lambda)){
                    mixing_probs <- optimParallel(par=probs, fn=qre_system, method="L-BFGS-B", lower=0, upper=1, lambda=lambda, payoffs_SIR=fallback_solved_values, problem_type=problem_type, choices=choice_vector, SIR_state=perceived_state, exog_parms=exog_parms)$par
                }
                if(!is.finite(lambda)){
                    mixing_probs <- c(1,0,0, 0,1,0, 0,0,1)
                }      
                qre.tm.vec[i] <- as.numeric(round(proc.time()[3] - qre.tm,2))
                message("Fallback QRE solved! Time taken: ", qre.tm.vec[i])

                # print(as.matrix(mixing_probs[1:3],nrow=1))
                # print(as.matrix(choice_vector,ncol=1))
                # as.numeric(choice_vector)

                l_S <- sum(mixing_probs[1:3]*choice_vector)
                l_I <- sum(mixing_probs[4:6]*choice_vector)
                l_R <- sum(mixing_probs[7:9]*choice_vector)

                # print(l_S)
                # message("yeetfeet")
            }
        }
        

        c_S = w*l_S
        c_I = w*phi*l_I
        c_R = w*l_R

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

    output_series <- cbind(output_series,qre_times=qre.tm.vec,mixing_prob_mat)

    if("eqm_concept"%in%names(parms)==TRUE) { 
        stopCluster(cl)
    }

    return(output_series)
}


# SIR_series_policy <- function(parms,choices_list,aggregates,policy_list) {
#     w = parms$w
#     phi = parms$phi
#     final.time = parms$final.time
#     pi_r = parms$pi_r
#     pi_d = parms$pi_d
#     rho_c = parms$rho_c
#     rho_l = parms$rho_l
#     rho_o = parms$rho_o
#     kappa_c = parms$kappa_c
#     kappa_l = parms$kappa_l

#     S_0 = aggregates$S
#     I_0 = aggregates$I
#     R_0 = aggregates$R
#     D_0 = 1 - (S_0 + I_0 + R_0)

#     policy_switch = policy_list$switch
#     policy_start = policy_list$start_date
#     policy_stop = policy_list$end_date
#     policy_labor_restriction = policy_list$labor_supply_level
#     policy_state_dependent = policy_list$state_dependent

#     state = data.frame(S=S_0,I=I_0,R=R_0)
#     types = c("S","I","R")
#     labor_supply = data.frame(S=NA,I=NA,R=NA)
#     lifetime_utility = data.frame(S=NA,I=NA,R=NA)

#     l_S0 = project_choices("S", state, choices_list$grid_list, choices_list$choices)
#     l_I0 = project_choices("I", state, choices_list$grid_list, choices_list$choices)
#     l_R0 = project_choices("R", state, choices_list$grid_list, choices_list$choices)

#     restricted_level = l_S0*policy_labor_restriction

#     output_series = data.frame(S = rep(S_0,length.out=final.time),
#                                 I = rep(I_0,length.out=final.time),
#                                 R = rep(R_0,length.out=final.time),
#                                 D = rep(D_0,length.out=final.time),
#                                 labor_S = rep(l_S0,length.out=final.time),
#                                 labor_I = rep(l_I0,length.out=final.time),
#                                 labor_R = rep(l_R0,length.out=final.time),
#                                 consumption_S = rep(w*l_S0,length.out=final.time),
#                                 consumption_I = rep(w*phi*l_I0,length.out=final.time),
#                                 consumption_R = rep(w*l_R0,length.out=final.time),
#                                 aggregate_labor_supply = rep((l_S0*S_0 + l_I0*I_0 + l_R0*R_0),length.out=final.time),
#                                 aggregate_consumption = rep((w*l_S0*S_0 + w*phi*l_I0*I_0 + w*l_R0*R_0),length.out=final.time)) 

#     for(i in 2:final.time) {
#         for(j in 1:length(types)) {
#             choices = choices_list$choices
#             agent_type = types[j]
#             labor_supply[j] = project_choices(agent_type, state, choices_list$grid_list, choices)
#             if((policy_switch==1)&(i>=policy_start)&(i<=policy_stop)) {
#                 labor_supply[j] = min(labor_supply[j],restricted_level)
#             }
#         }

#         l_S = labor_supply$S
#         l_I = labor_supply$I

#         c_S = w*l_S
#         c_I = w*phi*l_I

#         choices = data.frame(c_S = c_S, c_I = c_I, l_S = l_S, l_I = l_I)

#         state = SIR_update(parms,choices,state)
#         choice_vector = data.frame(
#             labor_S = labor_supply$S, 
#             labor_I = labor_supply$I, 
#             labor_R = labor_supply$R,
#             consumption_S = w*labor_supply$S,
#             consumption_I = w*phi*labor_supply$I,
#             consumption_R = w*labor_supply$R)
#         aggregate_labor_supply = labor_supply$S*state$S + labor_supply$I*state$I + labor_supply$R*state$R
#         aggregate_consumption = choice_vector$consumption_S*state$S + choice_vector$consumption_I*state$I + choice_vector$consumption_R*state$R

#         output_series[i,] = c(state, choice_vector, aggregate_labor_supply, aggregate_consumption)
#     }

#     Reff_choices = data.frame(l_S=output_series$labor_S, l_I=output_series$labor_I, l_R=output_series$labor_R, c_S=output_series$consumption_S, c_I=output_series$consumption_I, c_R=output_series$consumption_R)
#     Reff_SIR = data.frame(S=output_series$S, I=output_series$I, R=output_series$R, D=output_series$D)
#     Reff = compute_R0(parms, Reff_choices, Reff_SIR)

#     output_series = cbind(output_series, 
#                 weekly_new_cases = newly_infected(parms, Reff_choices, Reff_SIR),
#                 newly_recovered = pi_r*output_series$I,
#                 newly_dead = pi_d*output_series$I,
#                 newly_infected_per_I = newly_infected(parms, Reff_choices, Reff_SIR)/output_series$I,
#                 prob_infection_S = newly_infected(parms, Reff_choices, Reff_SIR)/output_series$S,
#                 infection_growth_rate = newly_infected(parms, Reff_choices, Reff_SIR)/output_series$I - pi_r - pi_d,
#                 Reff = Reff,
#                 consumption_contacts = rho_c*(output_series$consumption_S*output_series$consumption_I)^kappa_c,
#                 labor_contacts = rho_l*(output_series$labor_S*output_series$labor_I)^kappa_l,
#                 other_contacts = rho_o)
    
#     output_series = cbind(output_series, total_cases = cumsum(output_series$weekly_new_cases), aggregate_hours_deviation = 100*(output_series$aggregate_labor_supply/output_series$aggregate_labor_supply[1] - 1), aggregate_consumption_deviation = 100*(output_series$aggregate_consumption/output_series$aggregate_consumption[1] - 1))

#     return(output_series)
# }



# # function to generate important epi and econ series from solved policy functions and initial conditions for impulse response plots
# SIR_series_IRF <- function(parms,choices_list,aggregates,method="multilinear",grid_dfrm=NULL,impulse_parm=NULL,impulse_init=4,impulse_new=2,impulse_time=NULL, solved_values=NULL) {
#     w = parms$w
#     phi = parms$phi
#     final.time = parms$final.time
#     pi_r = parms$pi_r
#     pi_d = parms$pi_d
#     rho_c = impulse_init/(parms$Cbm)^2
#     rho_l = parms$rho_l
#     rho_o = parms$rho_o

#     # S_0 = aggregates$S
#     # I_0 = aggregates$I
#     # R_0 = aggregates$R
#     S_0 = aggregates[[1]]
#     I_0 = aggregates[[2]]
#     R_0 = aggregates[[3]]
#     D_0 = 1 - (S_0 + I_0 + R_0)

#     # state = c(S=S_0,I=I_0,R=R_0,benchmark_cons_contacts=impulse_init)
#     state = c(S=S_0,I=I_0,R=R_0)
#     types = c("S","I","R")
#     labor_supply = data.frame(S=NA,I=NA,R=NA)
#     lifetime_utility = data.frame(S=NA,I=NA,R=NA)

#     l_S = solved_values$labor_supply_S[which(solved_values$benchmark_cons_contacts==impulse_init)]
#     l_I = solved_values$labor_supply_I[which(solved_values$benchmark_cons_contacts==impulse_init)]
#     l_R = solved_values$labor_supply_R[which(solved_values$benchmark_cons_contacts==impulse_init)]

#     l_S_fn = ipol(l_S, grid=grid_list, method="multilinear")
#     l_I_fn = ipol(l_I, grid=grid_list, method="multilinear")
#     l_R_fn = ipol(l_R, grid=grid_list, method="multilinear")
#     l_S0 = l_S_fn(state)
#     l_I0 = l_I_fn(state)
#     l_R0 = l_R_fn(state)

#     output_series = data.frame(S = rep(S_0,length.out=final.time),
#                                 I = rep(I_0,length.out=final.time),
#                                 R = rep(R_0,length.out=final.time),
#                                 D = rep(D_0,length.out=final.time),
#                                 labor_S = rep(l_S0,length.out=final.time),
#                                 labor_I = rep(l_I0,length.out=final.time),
#                                 labor_R = rep(l_R0,length.out=final.time),
#                                 consumption_S = rep(w*l_S0,length.out=final.time),
#                                 consumption_I = rep(w*phi*l_I0,length.out=final.time),
#                                 consumption_R = rep(w*l_R0,length.out=final.time),
#                                 aggregate_labor_supply = rep((l_S0*S_0 + l_I0*I_0 + l_R0*R_0),length.out=final.time),
#                                 aggregate_consumption = rep((w*l_S0*S_0 + w*phi*l_I0*I_0 + w*l_R0*R_0),length.out=final.time)) 

#     print(l_S_fn(c(0.99,0.01,0)))
#     print(parms)

#     choices = choices_list$choices
#     for(i in 2:(impulse_time-1)) {
#         if(method=="multilinear") {
#             for(j in 1:length(types)) {
#                 # if(i>=impulse_time) state[4] = impulse_new
                
#                 agent_type = types[j]
#                 # print(state)
#                 state = as.numeric(state)
#                 if(agent_type=="S") labor_supply[j] = l_S_fn(state)
#                 if(agent_type=="I") labor_supply[j] = l_I_fn(state)
#                 if(agent_type=="R") labor_supply[j] = l_R_fn(state)
#                 # print(paste0("worked ",i))
#                 # print(labor_supply$S)
#             }
#             # print("Collect 200")
#         }

#         l_S = labor_supply$S
#         l_I = labor_supply$I

#         # print("Collect 400")

#         c_S = w*l_S
#         c_I = w*phi*l_I

#         choices = data.frame(c_S = c_S, c_I = c_I, l_S = l_S, l_I = l_I)

#         # state_old = data.frame(S = state[1], I = state[2], R = state[3], benchmark_cons_contacts = state[4])
#         # state = SIR_update(parms,choices,state_old)
#         # print(state)
#         state = data.frame(S = state[1], I = state[2], R = state[3])
#         state = SIR_update(parms,choices,state)
#         state = state[1:3]
#         # state[4] = state_old[4]
#         # print(state)
#         choice_vector = data.frame(
#             labor_S = labor_supply$S, 
#             labor_I = labor_supply$I, 
#             labor_R = labor_supply$R,
#             consumption_S = w*labor_supply$S,
#             consumption_I = w*phi*labor_supply$I,
#             consumption_R = w*labor_supply$R)
#         aggregate_labor_supply = labor_supply$S*state$S + labor_supply$I*state$I + labor_supply$R*state$R
#         aggregate_consumption = choice_vector$consumption_S*state$S + choice_vector$consumption_I*state$I + choice_vector$consumption_R*state$R

#         output_series[i,] = c(state, choice_vector, aggregate_labor_supply, aggregate_consumption)
#     }

#     l_S = solved_values$labor_supply_S[which(solved_values$benchmark_cons_contacts==impulse_new)]
#     l_I = solved_values$labor_supply_I[which(solved_values$benchmark_cons_contacts==impulse_new)]
#     l_R = solved_values$labor_supply_R[which(solved_values$benchmark_cons_contacts==impulse_new)]
#     l_S_fn = ipol(l_S, grid=grid_list, method="multilinear")
#     l_I_fn = ipol(l_I, grid=grid_list, method="multilinear")
#     l_R_fn = ipol(l_R, grid=grid_list, method="multilinear")
#     parms$rho_c = impulse_new/(parms$Cbm)^2

#     print(l_S_fn(c(0.99,0.01,0)))

#     for(i in impulse_time:final.time) {
#         if(method=="multilinear") {
#             for(j in 1:length(types)) {
                
#                 agent_type = types[j]
#                 state = as.numeric(state)
#                 if(agent_type=="S") labor_supply[j] = l_S_fn(state)
#                 if(agent_type=="I") labor_supply[j] = l_I_fn(state)
#                 if(agent_type=="R") labor_supply[j] = l_R_fn(state)
#             }
#         }

#         l_S = labor_supply$S
#         l_I = labor_supply$I

#         c_S = w*l_S
#         c_I = w*phi*l_I

#         choices = data.frame(c_S = c_S, c_I = c_I, l_S = l_S, l_I = l_I)

#         state = data.frame(S = state[1], I = state[2], R = state[3])
#         state = SIR_update(parms,choices,state)
#         state = state[1:3]
#         choice_vector = data.frame(
#             labor_S = labor_supply$S, 
#             labor_I = labor_supply$I, 
#             labor_R = labor_supply$R,
#             consumption_S = w*labor_supply$S,
#             consumption_I = w*phi*labor_supply$I,
#             consumption_R = w*labor_supply$R)
#         aggregate_labor_supply = labor_supply$S*state$S + labor_supply$I*state$I + labor_supply$R*state$R
#         aggregate_consumption = choice_vector$consumption_S*state$S + choice_vector$consumption_I*state$I + choice_vector$consumption_R*state$R

#         output_series[i,] = c(state, choice_vector, aggregate_labor_supply, aggregate_consumption)
#     }

#     Reff_choices = data.frame(l_S=output_series$labor_S, l_I=output_series$labor_I, l_R=output_series$labor_R, c_S=output_series$consumption_S, c_I=output_series$consumption_I, c_R=output_series$consumption_R)
#     Reff_SIR = data.frame(S=output_series$S, I=output_series$I, R=output_series$R, D=output_series$D)
#     Reff = compute_R0(parms, Reff_choices, Reff_SIR)

#     output_series = cbind(output_series, 
#                     weekly_new_cases = newly_infected(parms, Reff_choices, Reff_SIR),
#                     newly_recovered = pi_r*output_series$I,
#                     newly_dead = pi_d*output_series$I,
#                     newly_infected_per_I = newly_infected(parms, Reff_choices, Reff_SIR)/output_series$I,
#                     prob_infection_S = newly_infected(parms, Reff_choices, Reff_SIR)/output_series$S,
#                     infection_growth_rate = newly_infected(parms, Reff_choices, Reff_SIR)/output_series$I - pi_r - pi_d,
#                     Reff = Reff,
#                     consumption_contacts = rho_c*output_series$consumption_S*output_series$consumption_I,
#                     labor_contacts = rho_l*output_series$labor_S*output_series$labor_I,
#                     other_contacts = rho_o)
#     output_series = cbind(output_series, total_cases = cumsum(output_series$weekly_new_cases), aggregate_hours_deviation = 100*(output_series$aggregate_labor_supply/output_series$aggregate_labor_supply[1] - 1) , aggregate_consumption_deviation = 100*(output_series$aggregate_consumption/output_series$aggregate_consumption[1] - 1) )

#     return(output_series)
# }

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

# # function to calculate logit quantal response function with three player types. probs should be a 2kx1 vector and payoffs should be a kx3k matrix. k is the number of choices (3 -- l_S, l_I, l_R). Compute the payoffs assuming everyone's choices from the pre_value_fn. payoffs should be an object with named rows and columns.
lqrf <- function(lambda, payoffs, probs, action_type) {

    probs_1 <- matrix(probs[1:3]) # column vector
    probs_2 <- t(matrix(probs[4:6])) # row vector

    probs_mat <- probs_1%*%probs_2

    probs_num <- as.vector(probs_mat)
    payoffs_num <- as.matrix(payoffs)

    # print(sum(probs_num))

    expected_payoffs <- payoffs_num%*%probs_num
    noised_EP <- lambda*expected_payoffs
    exp_EU <- exp(noised_EP)

    evaluated_action <- which(rownames(payoffs)==action_type)
    # print(evaluated_action)

    # output <- exp(noised_EP[evaluated_action])/sum(exp(noised_EP)) # brute calc -- vulnerable to overflow. payoff differences are better here, less vulnerable

    payoff_diffs_wrt_evaluated_action <- as.numeric(noised_EP) - as.numeric(noised_EP[evaluated_action])
    den <- sum(exp(payoff_diffs_wrt_evaluated_action))


    output <- 1/den
    if(is.na(den==TRUE)) {output <- 0}

    # print(output)

    return(output)
}

# lambda <- 1000
# payoffs <- data.frame(I_playing_S = c(5,4), I_playing_I = c(3,3))
# rownames(payoffs) = c("S_playing_S", "S_playing_I")
# probs <- c(0.9,0.1)
# action_type <- "S_playing_S"
# lqrf(lambda, payoffs, probs, action_type)

# # solver to calculate a QRE given: an initial vector of choice probabilities "choice probs" (1x9), a noise parameter "lambda" (1x1), a matrix with the solved value functions for each type "payoffs_SIR" (ngridpts x 6, 3 for the states and 3 for the values), choices being mixed between "choices" (1x3), and current state of the system "SIR_state". Outputs a 1x9 vector of imbalances. When plugged into nleqslv it produces a 1x9 vector of mixing probabilities (1x3 for each type, across all three type-specific actions).
qre_system <- function(choice_probs, lambda, payoffs_SIR, problem_type="eqm", choices, SIR_state, exog_parms) {
    ### set up initial choice probabilities
    #### S
    S_choice_probs_S <- choice_probs[1]
    S_choice_probs_I <- choice_probs[2]
    S_choice_probs_R <- choice_probs[3]
    S_probs <- c(S_choice_probs_S, S_choice_probs_I, S_choice_probs_R)
    #### I
    I_choice_probs_S <- choice_probs[4]
    I_choice_probs_I <- choice_probs[5]
    I_choice_probs_R <- choice_probs[6]
    I_probs <- c(I_choice_probs_S, I_choice_probs_I, I_choice_probs_R)
    #### R
    R_choice_probs_S <- choice_probs[7]
    R_choice_probs_I <- choice_probs[8]
    R_choice_probs_R <- choice_probs[9]
    R_probs <- c(R_choice_probs_S, R_choice_probs_I, R_choice_probs_R)

    #### others' choice probs to go into the lqrfs. each should be a 6x1 vector which contains the 3x1 vectors for the other two types.
    SR_choice_probs <- c(S_choice_probs_S, S_choice_probs_I, S_choice_probs_R, R_choice_probs_S, R_choice_probs_I, R_choice_probs_R)
    IR_choice_probs <- c(I_choice_probs_S, I_choice_probs_I, I_choice_probs_R, R_choice_probs_S, R_choice_probs_I, R_choice_probs_R)
    SI_choice_probs <- c(S_choice_probs_S, S_choice_probs_I, S_choice_probs_R, I_choice_probs_S, I_choice_probs_I, I_choice_probs_R)

    ### construct mixed strategies using choice probabilities
    l_S <- as.matrix(S_probs)%*%as.matrix(choices)
    l_I <- as.matrix(I_probs)%*%as.matrix(choices)
    l_R <- as.matrix(R_probs)%*%as.matrix(choices)

    ## compute payoff matrix for each type
    ### In the decentralized problems, people use their own value functions to choose their mixtures. The LQRF here represents that they are uncertain over their state, and so their "signal" regarding which value function to use is contaminated by noise (since the test isn't a strong enough signal to overwhelm the noise). They have some private information which helps them determine their type, and the better the test is the more effectively they're able to combine it with their private information to determine their type and the corresponding appropriate action. What's interesting is that there's a test quality < 1 where the S types decide the evidence is strong enough and play their proper action.
    if(problem_type=="eqm"){ 
        payoffs_S <- compute_payoff_matrix(choices, payoffs=payoffs_SIR, exog_parms=exog_parms, SIR_state=SIR_state, grid_list=grid_list, agent_type="S") 
        payoffs_I <- compute_payoff_matrix(choices, payoffs=payoffs_SIR, exog_parms=exog_parms, SIR_state=SIR_state, grid_list=grid_list, agent_type="I") 
        payoffs_R <- compute_payoff_matrix(choices, payoffs=payoffs_SIR, exog_parms=exog_parms, SIR_state=SIR_state, grid_list=grid_list, agent_type="R") 
    }
    ### Under the coordinated problem, we've solved the incentives and people want to do whatever maximizes the SWF -- they just don't know what their type is, and so have trouble determining which action from them would maximize the SWF.
    if(problem_type=="plan"){ 
        payoffs_S <- compute_payoff_matrix(choices, payoffs=payoffs_SIR, exog_parms=exog_parms, SIR_state=SIR_state, grid_list=grid_list, agent_type="S", problem_type="plan") 
        payoffs_I <- compute_payoff_matrix(choices, payoffs=payoffs_SIR, exog_parms=exog_parms, SIR_state=SIR_state, grid_list=grid_list, agent_type="I", problem_type="plan") 
        payoffs_R <- compute_payoff_matrix(choices, payoffs=payoffs_SIR, exog_parms=exog_parms, SIR_state=SIR_state, grid_list=grid_list, agent_type="R", problem_type="plan") 
    }


    rownames(payoffs_S) <- c("S_playing_S", "S_playing_I", "S_playing_R")
    rownames(payoffs_I) <- c("I_playing_S", "I_playing_I", "I_playing_R")
    rownames(payoffs_R) <- c("R_playing_S", "R_playing_I", "R_playing_R")

    # print(payoffs_S)
    # print(payoffs_I)
    # print(payoffs_R)

    ### calculate lqrf choice probabilities
    # if(is.finite(lambda)){
        #### S
        S_choice_probs_S_update <- lqrf(lambda, payoffs_S, IR_choice_probs, "S_playing_S")
        S_choice_probs_I_update <- lqrf(lambda, payoffs_S, IR_choice_probs, "S_playing_I")
        S_choice_probs_R_update <- 1 - S_choice_probs_S_update - S_choice_probs_I_update
        #### I
        I_choice_probs_S_update <- lqrf(lambda, payoffs_I, SR_choice_probs, "I_playing_S")
        I_choice_probs_I_update <- lqrf(lambda, payoffs_I, SR_choice_probs, "I_playing_I")
        # print(I_choice_probs_S_update)
        # print(I_choice_probs_I_update)
        I_choice_probs_R_update <- 1 - I_choice_probs_S_update - I_choice_probs_I_update
        #### R
        R_choice_probs_S_update <- lqrf(lambda, payoffs_R, SI_choice_probs, "R_playing_S")
        R_choice_probs_I_update <- lqrf(lambda, payoffs_R, SI_choice_probs, "R_playing_I")
        R_choice_probs_R_update <- 1 - R_choice_probs_S_update - R_choice_probs_I_update
    # }

    # if(!is.finite(lambda)) {
    #     S_choice_probs_S_update <- 1
    #     S_choice_probs_I_update <- 0
    #     S_choice_probs_R_update <- 0

    #     I_choice_probs_S_update <- 0
    #     I_choice_probs_I_update <- 1
    #     I_choice_probs_R_update <- 0

    #     R_choice_probs_S_update <- 0
    #     R_choice_probs_I_update <- 0
    #     R_choice_probs_R_update <- 1
    # }

    # print(c(S_choice_probs_S_update, S_choice_probs_I_update, S_choice_probs_R_update,
        # I_choice_probs_S_update, I_choice_probs_I_update, I_choice_probs_R_update,
        # R_choice_probs_S_update, R_choice_probs_I_update, R_choice_probs_R_update))

    ### calculate squared imbalances
    #### S
    S_imbalance_S <- (S_choice_probs_S - S_choice_probs_S_update)^2 
    S_imbalance_I <- (S_choice_probs_I - S_choice_probs_I_update)^2
    S_imbalance_R <- (S_choice_probs_R - S_choice_probs_R_update)^2
    #### I
    I_imbalance_S <- (I_choice_probs_S - I_choice_probs_S_update)^2 
    I_imbalance_I <- (I_choice_probs_I - I_choice_probs_I_update)^2
    I_imbalance_R <- (I_choice_probs_R - I_choice_probs_R_update)^2
    #### R
    R_imbalance_S <- (R_choice_probs_S - R_choice_probs_S_update)^2 
    R_imbalance_I <- (R_choice_probs_I - R_choice_probs_I_update)^2
    R_imbalance_R <- (R_choice_probs_R - R_choice_probs_R_update)^2
    ### total
    total_imbalance <- S_imbalance_S + S_imbalance_I + S_imbalance_R + I_imbalance_S + I_imbalance_I + I_imbalance_R + R_imbalance_S + R_imbalance_I + R_imbalance_R
    # probability_penalty <- (1 - (S_choice_probs_S_update + S_choice_probs_I_update + S_choice_probs_R_update))^2 + (1 - (I_choice_probs_S_update + I_choice_probs_I_update + I_choice_probs_R_update))^2 + (1 - (R_choice_probs_S_update + R_choice_probs_I_update + R_choice_probs_R_update))^2 

    # print(total_imbalance)
    # print(probability_penalty)
    loss <- total_imbalance
    # loss <- total_imbalance + probability_penalty

    # return(c(S_imbalance_S, S_imbalance_I, S_imbalance_R, I_imbalance_S, I_imbalance_I, I_imbalance_R, R_imbalance_S, R_imbalance_I, R_imbalance_R)) # use with nleqslv
    return(loss) # use with optim
}

# lambda <- 1
# probs <- c(1,0,0, 0,1,0, 0,0,1)
# qre_system(choice_probs=probs, lambda=lambda, payoffs_SIR=eqm_vpfn, choices=data.frame(S = 7, I = 4, R = 8), SIR_state=SIR_init, exog_parms=exog_parms)

# # nleqslv(x=probs, fn=qre_system, lambda=lambda, payoffs_SIR=eqm_vpfn, choices=data.frame(S = 7, I = 4, R = 8), SIR_state=SIR_init, exog_parms=exog_parms)

# optim(par=probs, fn=qre_system, method="L-BFGS-B", lower=0, upper=1, lambda=100, payoffs_SIR=eqm_vpfn, choices=data.frame(S = 7.99, I = 7.72, R = 7.99), SIR_state=SIR_init, exog_parms=exog_parms)
# optim(par=probs, fn=qre_system, method="L-BFGS-B", lower=0, upper=1, lambda=100, payoffs_SIR=eqm_vpfn, choices=data.frame(S = eqm_vpfn$labor_supply_S[50], I = eqm_vpfn$labor_supply_I[50], R = eqm_vpfn$labor_supply_R[50]), SIR_state=data.frame(S=eqm_vpfn$S[50],I=eqm_vpfn$I[50],R=eqm_vpfn$R[50] ), exog_parms=exog_parms)

# Computes the payoff matrix for an agent from choosing one of three actions given all the possibilities for the other players. "own_choice" is a 1x1 of what the current type will choose, "others_choices" is a 1x2 vector of what the others will choose. "payoffs" should be a long nx4 of the own-type value function at each state (3 for the state, 1 for the value). "exog_parms" is the vector of exogenous parameters, "SIR_state" is a 1x3 of the current disease state, "agent_type" is a 1x1 character from {S,I,R}, "grid_list" is a list of the grid used for the payoffs, "contval_list" is a list of the value function. The output of this should be a 3x9 matrix.
# The choices going in should already be probability-weighted. No probability-weighting happens inside here. "choices" should be a dataframe with three named columns, S I R, and be 1x3.
compute_payoff_matrix <- function(choices, payoffs, exog_parms, SIR_state, grid_list, agent_type, problem_type="eqm") {

    contval_list = list(S=payoffs$lifetime_utility_S, I=payoffs$lifetime_utility_I, R=payoffs$lifetime_utility_R, SWF=payoffs$SWF)

    all_choices <- data.frame(S=NA,I=NA,R=NA)

    U_ijS <- matrix(NA,nrow=3,ncol=3)
    U_ijI <- matrix(NA,nrow=3,ncol=3)
    U_ijR <- matrix(NA,nrow=3,ncol=3)

    for(i in 1:ncol(choices)){
        for(j in 1:ncol(choices)){
            for(k in 1:ncol(choices)) {
                all_choices$S <- as.numeric(choices[i])
                all_choices$I <- as.numeric(choices[j])
                all_choices$R <- as.numeric(choices[k])
                if(problem_type=="eqm"){
                    choice <- as.numeric(all_choices[agent_type])
                    if(k==1) {U_ijS[i,j] <- pre_value_fn(choice, exog_parms, all_choices, SIR_state, agent_type, grid_list, contval_list)}
                    if(k==2) {U_ijI[i,j] <- pre_value_fn(choice, exog_parms, all_choices, SIR_state, agent_type, grid_list, contval_list)}
                    if(k==3) {U_ijR[i,j] <- pre_value_fn(choice, exog_parms, all_choices, SIR_state, agent_type, grid_list, contval_list)}
                }
                if(problem_type=="plan"){
                    if(k==1) {U_ijS[i,j] <- pre_value_fn_planner(all_choices, exog_parms, SIR_state,grid_list, contval_list)}
                    if(k==2) {U_ijI[i,j] <- pre_value_fn_planner(all_choices, exog_parms, SIR_state,grid_list, contval_list)}
                    if(k==3) {U_ijR[i,j] <- pre_value_fn_planner(all_choices, exog_parms, SIR_state,grid_list, contval_list)}
                }
                
            } 
        } 
    }

    if(agent_type=="S") { U_ijk <- as.matrix(cbind(U_ijS, U_ijI, U_ijR)) }
    if(agent_type=="I") { U_ijk <- as.matrix(cbind(t(U_ijS), t(U_ijI), t(U_ijR))) }
    if(agent_type=="R") { U_ijk <- as.matrix(rbind( as.vector(U_ijS), as.vector(U_ijI), as.vector(U_ijR) )) }

    return(U_ijk)
}

# choices <- data.frame(S = eqm_vpfn$labor_supply_S[50], I = eqm_vpfn$labor_supply_I[50], R = eqm_vpfn$labor_supply_R[50])
# choices <- data.frame(S = 7, I = 4, R = 8)
# SIR_state <- data.frame(S=0.49, I=0.01, R=0.5)

# compute_payoff_matrix(choices, payoffs=eqm_vpfn, exog_parms=exog_parms, SIR_state=SIR_state, grid_list=grid_list, agent_type="S") 

# probs <- c(0.8,0.1,0.1, 0.1,0.8,0.1, 0.1,0.1,0.8)
# payoffs_S <- data.frame(I_playing_S = c(5,4,3), I_playing_I = c(3,3,3), I_playing_R = c(3,3,3))
# rownames(payoffs_S) = c("S_playing_S", "S_playing_I", "S_playing_R")
# payoffs_I <- data.frame(S_playing_I = c(3,5,3), S_playing_S = c(4,5,3), S_playing_R = c(2,2,2))
# rownames(payoffs_I) = c("I_playing_S", "I_playing_I", "I_playing_R")
# payoffs_S
# payoffs_I

# nleqslv(probs, qre_system, payoffs_S = payoffs_S, payoffs_I = payoffs_I, lambda = 1000)

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

    c_S = as.numeric(w*l_S)
    c_I = as.numeric(w*phi*l_I)
    c_R = as.numeric(w*l_R)
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

    caseload_wt <- exog_parms$caseload_wt

    plannertype = exog_parms$plannertype

    objective_type <- exog_parms$objective_type

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

    if(exog_parms$caseload_type=="linear") {caseload_cost <- 1 - caseload_wt*I}
    if(exog_parms$caseload_type=="convex") {caseload_cost <- 1 - caseload_wt*exp(I)}

    value_S = total_utility(c_S,l_S,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha_U,Lbar)*caseload_cost + discount_factor*((1-tau_jt)*continuation_value(SIR_, grid_list, contval_list$S) + tau_jt*continuation_value(SIR_, grid_list, contval_list$I))

    value_I = total_utility(c_I,l_I,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha_U,Lbar)*caseload_cost + discount_factor*((1-pi_r-pi_d)*continuation_value(SIR_, grid_list, contval_list$I) + pi_r*continuation_value(SIR_, grid_list, contval_list$R)) + pi_d*U_D

    value_R = total_utility(c_R,l_R,C,L,gamma_c,gamma_l,rho_c,rho_l,risk_aversion,s,alpha_U,Lbar) + discount_factor*(continuation_value(SIR_, grid_list, contval_list$R))

    # if(plannertype=="unweighted") {value = value_S + value_I + value_R}
    # if(plannertype=="weighted") {value = S*value_S + I*value_I + R*value_R}
    if(objective_type=="welfare") {value = S*value_S + I*value_I + R*value_R}
    if(objective_type=="cases") { value = SIR_$I }

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
generate_time_series <- function(solved_values, interpolator="multilinear", grid_dfrm=NULL, infolag_list=NULL, fallback_solved_values=NULL, policy_list=NULL, ...) {
    others_choices = data.frame(l_S = solved_values$labor_supply_S, l_I = solved_values$labor_supply_I, l_R = solved_values$labor_supply_R, simplex_check=simplex_check)
    
    others_utilities = data.frame(u_S = solved_values$lifetime_utility_S, u_I = solved_values$lifetime_utility_I, u_R = solved_values$lifetime_utility_R)

    choices_list = list(grid_list = grid_list, choices = others_choices, utilities = others_utilities)

    fallback_choices_list <- NULL
    if(is.null(fallback_solved_values)==FALSE) {
        # message("yeetbeets")
        fallback_others_choices = data.frame(l_S = fallback_solved_values$labor_supply_S, l_I = fallback_solved_values$labor_supply_I, l_R = fallback_solved_values$labor_supply_R, simplex_check=simplex_check)
        fallback_others_utilities = data.frame(u_S = fallback_solved_values$lifetime_utility_S, u_I = fallback_solved_values$lifetime_utility_I, u_R = fallback_solved_values$lifetime_utility_R)
        fallback_choices_list <- list(grid_list = grid_list, choices = fallback_others_choices, utilities = fallback_others_utilities)
    }
    projection = SIR_series(exog_parms, choices_list, SIR_init, method=interpolator, grid_dfrm=grid_dfrm, infolag_list=infolag_list, fallback_choices_list=fallback_choices_list, policy_list=policy_list, grid_list=grid_list, solved_values=solved_values, fallback_solved_values=fallback_solved_values)

    results = cbind(time=c(1:final.time),projection)

    return(results)
}

# # Generates a simple time series from an initial point, parameterization, and interpolator. Should output SIR, labor supply, consumption, and contact series.
# simple_SIR_series <- function(parms,initial_condition,interpolators,sim_time) {
#     w = parms$w
#     phi = parms$phi
#     final.time = parms$final.time
#     pi_r = parms$pi_r
#     pi_d = parms$pi_d
#     rho_c = parms$rho_c
#     rho_l = parms$rho_l
#     rho_o = parms$rho_o

#     S_0 = initial_condition[1]
#     I_0 = initial_condition[2]
#     R_0 = initial_condition[3]
#     D_0 = 1 - (S_0 + I_0 + R_0)

#     # state = c(S=S_0,I=I_0,R=R_0,benchmark_cons_contacts=impulse_init)
#     state = data.frame(S=S_0,I=I_0,R=R_0,D=D_0) #for contact() and compute_R0()
#     state_small = c(S=state[[1]],I=state[[2]],R=state[[3]]) #for interpolators
#     types = c("S","I","R")
#     labor_supply = data.frame(S=NA,I=NA,R=NA)
#     # lifetime_utility = data.frame(S=NA,I=NA,R=NA)

#     l_S = interpolators[["ls_fn"]]
#     l_I = interpolators[["li_fn"]]
#     l_R = interpolators[["lr_fn"]]

#     l_S0 = l_S(state_small)
#     l_I0 = l_I(state_small)
#     l_R0 = l_R(state_small)

#     choices_0 = data.frame(c_S = l_S(state_small)*w, c_I = l_I(state_small)*w*phi, l_S = l_S(state_small), l_I = l_I(state_small))

#     initial_contacts = contacts(parms, choices_0)
#     R0 = compute_R0(parms, choices_0, state)

#     output_series = matrix(-1,nrow=sim_time, ncol=14)
#     output_series[1,] = c(S_0,I_0,R_0,D_0,l_S0,l_I0,l_R0,w*l_S0,w*phi*l_I0,w*l_R0,(l_S0*S_0 + l_I0*I_0 + l_R0*R_0),(w*l_S0*S_0 + w*phi*l_I0*I_0 + w*l_R0*R_0),initial_contacts,R0)

#     # output_series = matrix(rep(S_0,length.out=sim_time),
#     #                             rep(I_0,length.out=sim_time),
#     #                             rep(R_0,length.out=sim_time),
#     #                             rep(D_0,length.out=sim_time),
#     #                             rep(l_S0,length.out=sim_time),
#     #                             rep(l_I0,length.out=sim_time),
#     #                             rep(l_R0,length.out=sim_time),
#     #                             rep(w*l_S0,length.out=sim_time),
#     #                             rep(w*phi*l_I0,length.out=sim_time),
#     #                             rep(w*l_R0,length.out=sim_time),
#     #                             rep((l_S0*S_0 + l_I0*I_0 + l_R0*R_0),length.out=sim_time),
#     #                             rep((w*l_S0*S_0 + w*phi*l_I0*I_0 + w*l_R0*R_0),length.out=sim_time),
#     #                             rep(initial_contacts,length.out=sim_time),
#     #                             rep(R0,length.out=sim_time), byrow=TRUE) 
#     output_series <- as.data.frame(output_series)
#     colnames(output_series) = c("S","I","R","D","labor_S","labor_I","labor_R","consumption_S", "consumption_I", "consumption_R", "aggregate_labor_supply", "aggregate_consumption", "contacts", "Reff")

#     state = SIR_update(parms, choices_0, state)

#     for(i in 2:sim_time) {
#         state_small = c(S=state[[1]],I=state[[2]],R=state[[3]])
#         labor_supply$S = l_S(state_small)
#         labor_supply$I = l_I(state_small)
#         labor_supply$R = l_R(state_small)

#         choices = data.frame(c_S = labor_supply$S*w, c_I = w*phi*labor_supply$I, l_S = labor_supply$S, l_I = labor_supply$I)
#         current_contacts = contacts(parms, choices)
#         aggregate_labor_supply = labor_supply$S*state_small[[1]] + labor_supply$I*state_small[[2]] + labor_supply$R*state_small[[3]]
#         aggregate_consumption = w*labor_supply$S*state_small[[1]] + w*phi*labor_supply$I*state_small[[2]] + w*labor_supply$R*state_small[[3]]

#         new_row = c(S = state$S,
#                                 I = state_small[[1]],
#                                 R = state_small[[2]],
#                                 D = state_small[[3]],
#                                 labor_S = labor_supply$S,
#                                 labor_I = labor_supply$I,
#                                 labor_R = labor_supply$R,
#                                 consumption_S = labor_supply$S*w,
#                                 consumption_I = labor_supply$I*w*phi,
#                                 consumption_R = labor_supply$R*w,
#                                 aggregate_labor_supply = aggregate_labor_supply,
#                                 aggregate_consumption = aggregate_consumption,
#                                 contacts = current_contacts,
#                                 Reff = compute_R0(parms,choices,state))

#         output_series[i,] = new_row

#         state = SIR_update(parms, choices, state)
#     }

#     return(output_series)
# }

# # Generates impulse response functions
# generate_IRF <- function(SIR_init, impulse_1_list, impulse_2_list, grid_list=NULL, impulse_time, final_time) {
    
#     exog_parms_1 = impulse_1_list$exog_parms
#     series_1 = simple_SIR_series(exog_parms_1, SIR_init, impulse_1_list, (impulse_time-1))

#     print(head(series_1))

#     new_init_cond = series_1[(impulse_time-1), 1:4]

#     exog_parms_2 = impulse_2_list$exog_parms
#     series_2 = simple_SIR_series(exog_parms_2, SIR_init, impulse_2_list, (final_time - impulse_time + 1))

#     series = rbind(series_1,series_2)
#     results = cbind(time=seq(1:final_time),series)

#     return(results)
# }

# # wrapper for SIR_series. generates a time series of SIR outcomes from a dataframe of policy functions
# generate_time_series_policy <- function(solved_values, policy_list, ...) {
#     others_choices = data.frame(l_S = solved_values$labor_supply_S, l_I = solved_values$labor_supply_I, l_R = solved_values$labor_supply_R)

#     others_utilities = data.frame(u_S = solved_values$lifetime_utility_S, l_I = solved_values$lifetime_utility_I, l_R = solved_values$lifetime_utility_R)

#     choices_list = list(grid_list = grid_list, choices = others_choices, utilities = others_utilities)
#     projection = SIR_series_policy(exog_parms, choices_list, SIR_init, policy_list)

#     results = cbind(time=c(1:final.time),projection)

#     return(results)
# }


# runs the dynamic programming solver loop
dp_solver <- function(...) {

    ## initialize solver deltas and set tolerances
    damping_factor = 1

    VFI_count = 1

    l_lowerbound = 1e-9
    l_upperbound = Lbar*0.99
    l_upperbound_uniroot = Lbar*0.99
    epsilon_VFI = abs(0.0075*mean(as.numeric(contval_list$S)))
    if(problemtype=="planner"){epsilon_VFI = abs(0.0001*mean(as.numeric(contval_list$SWF)))} # Too low a tolerance leaves artifacts in the S-type policy function.
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
                        # message("Finished consistency iteration ", iteration_consistency, ". Delta is ", round(delta_consistency,5), ", new labor supplies are (l_S, l_I, l_R) = (", round(l_S_new,3), ", ", round(l_I_new,3), ", ", round(l_R_new,3), "). \n")
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
            if( (delta_VFI_old - delta_VFI < min(1e-04,epsilon_VFI) ) & (delta_VFI_old - delta_VFI >= 0) & (VFI_count > 100) ) { delta_VFI = epsilon_VFI/2}
        }
        if(problemtype=="planner") {
            new_planval_mat = old_planval_mat

            new_planval_mat[simplex_check_1] = c(as.numeric(VFI_results[,7]))

            delta_VFI_old = delta_VFI
            delta_VFI = max(abs(new_planval_mat - old_planval_mat))
            if( (delta_VFI_old - delta_VFI < min(1e-04,epsilon_VFI) ) & (delta_VFI_old - delta_VFI >= 0) & (VFI_count > 100) ) { delta_VFI = epsilon_VFI/2}
            if( VFI_count > 120 ) { delta_VFI = epsilon_VFI/2}
            #if( (VFI_count > 2) ) { delta_VFI = 1e-04}
        }
        message("Finished iteration ", VFI_count, ". Delta is ", round(delta_VFI,10) ,". Change in delta is ", (delta_VFI_old - delta_VFI) ,". Total time: ", round(VFI_count_time,3) ," seconds. Average core-time per state: ", ncores*round(VFI_count_time/nrow(grid_dfrm),3), " seconds. Average value function levels: (S,I,R) = (", round(mean(VFI_results[,4]),1),",", round(mean(VFI_results[,5]),1), ",", round(mean(VFI_results[,6]),1),"), SWF = ", round(mean(VFI_results[,7]),1),".") 

        intermediate_results <- as.data.frame(cbind(grid_dfrm,VFI_results))
        colnames(intermediate_results) <- c("S","I","R","labor_supply_S","labor_supply_I","labor_supply_R","lifetime_utility_S","lifetime_utility_I","lifetime_utility_R","SWF")

        if(precision>=0) {fwrite(intermediate_results, file=paste0("../../Results/value_policy_functions/value_policy_functions__",scenario_label,".csv"))}

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

### Generates rows for scenarios in sensitivity figure
sensitivity_rows <- function(data, label) {

    Mix<-c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99") 
    label <- as.character(label)

    dynamics_plot_base <- ggplot(data=data, aes(x=time, group=as.character(type), color=as.character(type)))

    infection <- dynamics_plot_base + 
              geom_line(aes(y=I*100), size=1) +
              theme_bw() + ggtitle(paste0("Infecteds under ", label)) +
              labs(x="Day", y="Proportion of population (%)", color="Policy type") +
              # geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
              # geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
              scale_color_viridis(discrete=TRUE) +
              guides(linetype=FALSE, size=FALSE, color=FALSE)
    infection

    recession <- dynamics_plot_base + 
              geom_line(aes(y=aggregate_consumption_deviation), size=1) +
              theme_bw() + ggtitle(paste0("Recession under ", label)) +
              labs(x="Day", y="GDP deviation (%)", color="Policy type") +
              # geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
              # geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
              scale_color_viridis(discrete=TRUE) +
              guides(linetype=FALSE, size=FALSE)
    recession

    big_dfrm <- data %>% 
                    mutate(weighted_labor_S = S*labor_S) %>%
                    mutate(weighted_labor_I = I*labor_I) %>%
                    mutate(weighted_labor_R = R*labor_R) %>% 
                    mutate(total_contacts = consumption_contacts + labor_contacts + other_contacts) %>%
                    mutate(type = recode(type, eqm = "decentralized", plan = "coordinated")) %>%
                    mutate(SI_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_I, l_S = labor_S, l_I = labor_I) )) %>%
                    mutate(SR_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_R, l_S = labor_S, l_I = labor_R) )) %>%
                    mutate(RI_contacts = contacts(parms, choices = data.frame(c_S = consumption_R, c_I = consumption_I, l_S = labor_R, l_I = labor_I) )) %>%
                    mutate(SS_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_S, l_S = labor_S, l_I = labor_S) )) %>%
                    mutate(RR_contacts = contacts(parms, choices = data.frame(c_S = consumption_R, c_I = consumption_R, l_S = labor_R, l_I = labor_R) )) %>%
                    mutate(prob_contact_I = (SI_contacts + RI_contacts)/(SI_contacts + RI_contacts + SR_contacts + SS_contacts + RR_contacts)) %>%
                    mutate(prob_contact_I_weighted = (SI_contacts*S*I + RI_contacts*R*I)/(SI_contacts*S*I + RI_contacts*R*I + SR_contacts*S*R + SS_contacts*S*S + RR_contacts*R*R)) %>%
                    mutate(daily_new_cases_from_consumption = parms$tau*parms$rho_c*consumption_S*consumption_I*S*I) %>%
                    mutate(daily_new_cases_from_labor = parms$tau*parms$rho_l*labor_S*labor_I*S*I) %>%
                    mutate(daily_new_unavoidable_cases = parms$tau*parms$rho_o*S*I)

all_base <- ggplot(data=big_dfrm, aes(x=time))

weighted_labor_supplies = all_base +
            geom_line(aes(y = S*labor_S, group=type, linetype=type), size=1, color = "darkgreen") +
            geom_line(aes(y = I*labor_I, group=type, linetype=type), size=1, color = "firebrick4") +
            geom_line(aes(y = R*labor_R, group=type, linetype=type), size=1, color = "dodgerblue4") +
            theme_bw() + ggtitle("Aggregate labor supplies\n(green = S, red = I, blue = R)") +
            xlab("Day") + ylab("Normalized person-hours")

total_contacts = all_base + 
            geom_line(aes(y=total_contacts, group=type, linetype=type), size=1) +
            # geom_hline(yintercept = 12, linetype="dashed", color="darkgray") +
            # geom_hline(yintercept = 5, linetype="dashed", color="darkgray") +
            theme_bw() + ggtitle("Total S-I contacts") +
            xlab("Day") + ylab("Daily contacts") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

cases_by_site = all_base + 
            geom_line(aes(y=daily_new_cases_from_consumption, group=type, linetype=type), size=1, color="black") +
            geom_line(aes(y=daily_new_cases_from_labor, group=type, linetype=type), size=1, color="orange") +
            geom_line(aes(y=daily_new_unavoidable_cases, group=type, linetype=type), size=1, color="purple") +
            theme_bw() + ggtitle("Total cases by activity type\n(black=consumption, orange=labor, purple=other)") +
            xlab("Day") + ylab("Cases") #+ 
            #scale_colour_manual(values=Mix)

contacts_by_site = all_base + 
            geom_line(aes(y=consumption_contacts, group=type, linetype=type), size=1, color="black") +
            geom_line(aes(y=labor_contacts, group=type, linetype=type), size=1, color="orange") +
            geom_line(aes(y=other_contacts, group=type, linetype=type), size=1, color="purple") +
            theme_bw() + ggtitle("Average contacts by activity type\n(black=consumption, orange=labor, purple=other)") +
            xlab("Day") + ylab("Average contacts") #+ 
            #scale_colour_manual(values=Mix)

prob_contact_I_wt = all_base + 
            geom_line(aes(y=prob_contact_I_weighted, group=type, linetype=type), size=1) +
            theme_bw() + ggtitle("Probability random S or R contacts any I") +
            xlab("Day") + ylab("Probability") +
            scale_colour_manual(values=Mix) + guides(linetype=FALSE)


    scenarios_wide <- pivot_wider(big_dfrm[,c("time","type","aggregate_consumption","prob_contact_I_weighted","I","R","D","infolag")], id_cols=c("time","type","infolag"), names_from=type, values_from=c(aggregate_consumption, prob_contact_I_weighted, I, R, D), names_sep="_") %>%
                    mutate(total_implied_savings = aggregate_consumption_coordinated - aggregate_consumption_decentralized) %>%
                    mutate(PV_decentralized_losses = (58000/365 - aggregate_consumption_decentralized)*discount_factor^time) %>%
                    mutate(PV_coordinated_losses = (58000/365 -aggregate_consumption_coordinated)*discount_factor^time) %>%
                    mutate(PV_lockdown_losses = (58000/365 -aggregate_consumption_lockdown)*discount_factor^time) %>%
                    mutate(PV_total_implied_savings = (aggregate_consumption_coordinated - aggregate_consumption_decentralized)*discount_factor^time) %>%
                    mutate(total_averted_contacts_with_Is = (prob_contact_I_weighted_decentralized - prob_contact_I_weighted_coordinated) ) %>%
                    mutate(implied_savings_per_averted_I_contact = total_implied_savings/total_averted_contacts_with_Is)

    names(scenarios_wide)

    losses_and_deaths <- scenarios_wide %>% group_by(infolag) %>%
        summarise(decentralized_losses = round(sum(PV_decentralized_losses),2),
              decentralized_casesp100k = round(max(I_decentralized + R_decentralized + D_decentralized)*100000,0),
              coordinated_losses = round(sum(PV_coordinated_losses),2),
              coordinated_casesp100k = round(max(I_coordinated + R_coordinated + D_coordinated)*100000,0),
              lockdown_losses = round(sum(PV_lockdown_losses),2),
              lockdown_casesp100k = round(max(I_lockdown + R_lockdown + D_lockdown)*100000,0)
              )

    bar_data <- pivot_longer(losses_and_deaths, cols=c("decentralized_losses", "decentralized_casesp100k", "coordinated_losses", "coordinated_casesp100k", "lockdown_losses", "lockdown_casesp100k"), names_to = c("type","measure"), names_sep="_", values_to = "values")

    print(bar_data, nrow=100)

    bar_summary <- ggplot(bar_data, aes(fill=type, x=reorder(type, -values), y=values)) +
        geom_bar(position="dodge", stat="identity") +
        ggtitle(label) + 
        # facet_wrap(~'measure + infolag', scales = "free_y", labeller = as_labeller(c(casesp100k = "Total cases per 100,000",losses="Individual loss ($/person)")), strip.position = "left") +
        facet_wrap(~measure, labeller = as_labeller(c(casesp100k = "Total cases per 100,000",losses="Individual loss ($/person)")), strip.position = "left") +
       geom_text(aes(label=round(values,3)), position=position_dodge(width=0.9), vjust=-0.25) +
        #scale_fill_manual(values=Mix) +
        ylab(NULL) + xlab("") + labs(fill = "Policy type") +
        theme_classic() +
        theme(strip.background = element_blank(),
             strip.placement = "outside") +
        guides(fill=FALSE)

    output_list <- list(summary=bar_summary, infection=infection, recession=recession, weighted_labor_supplies=weighted_labor_supplies, cases_by_site=cases_by_site, contacts_by_site=contacts_by_site, prob_contact_I_wt=prob_contact_I_wt, total_contacts=total_contacts)

    return(output_list)
}

### Dynamics figures for information frictions cases
generate_dynamics_figures <- function(big_eqm_dfrm, big_opt_dfrm, label) {

      label <- as.character(label)

      eqm_plot_base <- ggplot(data=big_eqm_dfrm, aes(x=time, group=as.character(infolag), color=as.character(infolag)))
      opt_plot_base <- ggplot(data=big_opt_dfrm, aes(x=time, group=as.character(infolag), color=as.character(infolag)))

      infection_lower <- min(c(big_eqm_dfrm$I,big_opt_dfrm$I))*100
      infection_upper <- max(c(big_eqm_dfrm$I,big_opt_dfrm$I))*100

      supplies_lower <- min(c(big_eqm_dfrm$labor_S, big_eqm_dfrm$labor_I, big_eqm_dfrm$labor_R, big_opt_dfrm$labor_S, big_opt_dfrm$labor_I, big_opt_dfrm$labor_R))
      supplies_upper <- max(c(big_eqm_dfrm$labor_S, big_eqm_dfrm$labor_I, big_eqm_dfrm$labor_R, big_opt_dfrm$labor_S, big_opt_dfrm$labor_I, big_opt_dfrm$labor_R))

      eqm_infection <- eqm_plot_base + 
                  geom_line(aes(y=I*100), size=1) +
                  theme_bw() + ggtitle("Eqm: Infecteds") +
                  labs(x="Day", y="Proportion of population (%)", color=label) +
                  ylim(c(0,infection_upper)) + 
                  # geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
                  # geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
                  scale_color_viridis(discrete=TRUE) +
                  guides(linetype=FALSE, size=FALSE)
      eqm_infection

      eqm_labor_supplies <- eqm_plot_base + 
                  geom_line(aes(y=labor_S), size=1) +
                  geom_line(aes(y=labor_I), size=1.5, linetype = "dashed") +
                  ylim(c(supplies_lower,supplies_upper)) + 
                  theme_bw() + ggtitle("Eqm: S (solid) & I (dashed) labor supplies") +
                  labs(x="Day", y="Hours worked", color=label) +
                  # geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
                  # geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
                  scale_color_viridis(discrete=TRUE) +
                  guides(linetype=FALSE, size=FALSE)
      eqm_labor_supplies

      opt_infection <- opt_plot_base + 
                  geom_line(aes(y=I*100), size=1) +
                  theme_bw() + ggtitle("Plan: Infecteds") +
                  labs(x="Day", y="Proportion of population (%)", color=label) +
                  ylim(c(0,infection_upper)) + 
                  # geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
                  # geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
                  scale_color_viridis(discrete=TRUE) +
                  guides(linetype=FALSE, size=FALSE)
      opt_infection

      opt_labor_supplies <- opt_plot_base + 
                  geom_line(aes(y=labor_S), size=1) +
                  geom_line(aes(y=labor_I), size=1.5, linetype = "dashed") +
                  ylim(c(supplies_lower,supplies_upper)) + 
                  theme_bw() + ggtitle("Plan: S (solid) & I (dashed) labor supplies") +
                  labs(x="Day", y="Hours worked", color=label) +
                  # geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
                  # geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
                  scale_color_viridis(discrete=TRUE) +
                  guides(linetype=FALSE, size=FALSE)
      opt_labor_supplies

      lag_dynamics <- (eqm_infection / eqm_labor_supplies) | (opt_infection / opt_labor_supplies)

      return(lag_dynamics)
}

### Dynamics figures for information frictions cases -- single scenario
generate_dynamics_figures_small <- function(dfrm, label) {

      label <- as.character(label)

      plot_base <- ggplot(data=dfrm, aes(x=time, group=as.character(infolag), color=as.character(infolag)))

      infection <- plot_base + 
                  geom_line(aes(y=I*100), size=1) +
                  theme_bw() + ggtitle("Eqm: Infecteds") +
                  labs(x="Day", y="Proportion of population (%)", color=label) +
                  # ylim(c(0,infection_upper)) + 
                  # geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
                  # geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
                  scale_color_viridis(discrete=TRUE) +
                  guides(linetype=FALSE, size=FALSE, group=FALSE, color=FALSE)
      infection

      labor_supplies <- plot_base + 
                  geom_line(aes(y=labor_S), size=1) +
                  geom_line(aes(y=labor_I), size=1.5, linetype = "dashed") +
                  # ylim(c(supplies_lower,supplies_upper)) + 
                  theme_bw() + ggtitle("Eqm: S (solid) & I (dashed) labor supplies") +
                  labs(x="Day", y="Hours worked", color=label) +
                  # geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
                  # geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
                  scale_color_viridis(discrete=TRUE) +
                  guides(linetype=FALSE, size=FALSE, group=FALSE, color=FALSE)
      labor_supplies

      recession <- plot_base + 
                  geom_line(aes(y=aggregate_consumption_deviation), size=1) +
                  # ylim(c(supplies_lower,supplies_upper)) + 
                  theme_bw() + ggtitle("Recession") +
                  labs(x="Day", y="% deviation from initial level", color=label) +
                  # geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
                  # geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
                  scale_color_viridis(discrete=TRUE) +
                  guides(linetype=FALSE, size=FALSE, group=FALSE, color=FALSE)
      recession

      plotlist <- list(infection = infection, recession = recession, labor_supplies = labor_supplies)

      return(plotlist)
}

### Dynamics figures for information frictions cases -- combine three scenarios into a single panel
generate_dynamics_figures_small_composite <- function(dfrm, label, end_time=250, type="default") {

      label <- as.character(label)

      dfrm <- dfrm %>% filter(time <= end_time)

      write.csv(dfrm, file=paste0(label,"_dynamics.csv"))

      if(type=="default"){
            dfrm$type <- recode(dfrm$type, eqm = "\nVoluntary\nisolation\n", plan = "\nTargeted\nisolation\n", lockdown = "Lockdown")
            cols <- c("\nTargeted\nisolation\n" = viridis(3)[3], "\nVoluntary\nisolation\n" = viridis(3)[1],  "\nLockdown\n" = viridis(3)[2])
      }

      if(type=="eqm_plan"){
            dfrm$type <- recode(dfrm$type, eqm = "\nVoluntary\nisolation\n", plan = "\nTargeted\nisolation\n", lockdown_mech = "ld.casemin", lockdown_beh = "ld.decen", nocontrol = "\nNo control\n")
            cols <- c("\nTargeted\nisolation\n" = viridis(3)[3], "\nVoluntary\nisolation\n" = viridis(3)[1],  "\nNo control\n" = viridis(3)[2])
      }

      if(type=="case_penalty"){
            dfrm$type <- recode(dfrm$type, 
                                no.control = "\nNo control\n", 
                                penalty.0 = "\nTargeted\nisolation\n",
                                penalty.1000 = "\nTargeted\nisolation\nw/ case penalty\n")
            cols <- c("\nTargeted\nisolation\n" = viridis(3)[3], "\nTargeted\nisolation\nw/ case penalty\n" = viridis(3)[1],  "\nNo control\n" = viridis(3)[2])
      }

      plot_base <- ggplot(data=dfrm, aes(x=time, group=as.character(type), color=as.character(type)))

      infection <- plot_base + 
                    geom_line(aes(y=I*100), size=1) +
                    # geom_jitter(aes(y=I*100), size=1.2, width=1.2) +
                    theme_bw() + ggtitle(paste0("Infecteds ",label)) +
                    labs(x="Day", y="% currently\ninfected", color="Control type") +
                    ylim(c(0,max(dfrm$I)*100)) + 
                    # geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
                    # geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
                    scale_color_viridis(discrete=TRUE) +
                    theme(text = element_text(size=24)) +
                    guides(linetype=FALSE, size=FALSE, group=FALSE, color=FALSE)
      infection

      recession <- plot_base + 
                    geom_line(aes(y=aggregate_consumption_deviation), size=1) +
                    ylim(c(-100,0)) + 
                    theme_bw() + ggtitle(paste0("Recession ",label)) +
                    labs(x="Day", y="% deviation", color="Control type") +
                    # geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
                    # geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
                    scale_color_viridis(discrete=TRUE) +
                    theme(text = element_text(size=24)) +
                    guides(linetype=FALSE, size=FALSE, group=FALSE, color=FALSE)
      recession

      reff <- plot_base + 
                    geom_line(aes(y=Reff), size=1) +
                    theme_bw() + ggtitle(paste0("Reproductive number ",label)) +
                    labs(x="Day", y="R_eff", color="Control type") +
                    # geom_vline(xintercept = as.numeric(best_params["start_date",]), linetype="dashed") + 
                    # geom_vline(xintercept = as.numeric(best_params["end_date",]), linetype="dashed") + 
                    scale_color_viridis(discrete=TRUE) +
                    theme(text = element_text(size=24)) +
                    guides(linetype=FALSE, size=FALSE, group=FALSE)
      reff


    if(length(unique(dfrm$type))==3){ 
        infection <- infection + scale_color_manual(values=cols)  
        recession <- recession + scale_color_manual(values=cols)  
        reff <- reff + scale_color_manual(values=cols)  
    }

      plotlist <- list(infection = infection, recession = recession, reff = reff)

      return(plotlist)
}

# Function to generate recession metrics for information frictions scenarios
generate_recession_metrics <- function(big_eqm_dfrm, big_opt_dfrm, labels) {

      label <- as.character(labels[1])

      case1_label <- as.character(labels[2])
      case2_label <- as.character(labels[3])
      case3_label <- as.character(labels[4])

      big_dfrm <- rbind(big_eqm_dfrm, big_opt_dfrm)

      big_dfrm <- big_dfrm %>% 
                        mutate(weighted_labor_S = S*labor_S) %>%
                        mutate(weighted_labor_I = I*labor_I) %>%
                        mutate(weighted_labor_R = R*labor_R) %>% 
                        mutate(total_contacts = consumption_contacts + labor_contacts + other_contacts) %>%
                        mutate(type = recode(type, eqm = "decentralized", plan = "coordinated")) %>%
                        mutate(SI_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_I, l_S = labor_S, l_I = labor_I) )) %>%
                        mutate(SR_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_R, l_S = labor_S, l_I = labor_R) )) %>%
                        mutate(RI_contacts = contacts(parms, choices = data.frame(c_S = consumption_R, c_I = consumption_I, l_S = labor_R, l_I = labor_I) )) %>%
                        mutate(SS_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_S, l_S = labor_S, l_I = labor_S) )) %>%
                        mutate(RR_contacts = contacts(parms, choices = data.frame(c_S = consumption_R, c_I = consumption_R, l_S = labor_R, l_I = labor_R) )) %>%
                        mutate(prob_contact_I = (SI_contacts + RI_contacts)/(SI_contacts + RI_contacts + SR_contacts + SS_contacts + RR_contacts)) %>%
                        mutate(prob_contact_I_weighted = (SI_contacts*S*I + RI_contacts*R*I)/(SI_contacts*S*I + RI_contacts*R*I + SR_contacts*S*R + SS_contacts*S*S + RR_contacts*R*R)) %>%
                        mutate(daily_new_cases_from_consumption = parms$tau*parms$rho_c*consumption_S*consumption_I*S*I) %>%
                        mutate(daily_new_cases_from_labor = parms$tau*parms$rho_l*labor_S*labor_I*S*I) %>%
                        mutate(daily_new_unavoidable_cases = parms$tau*parms$rho_o*S*I)

      scenarios_wide <- pivot_wider(big_dfrm[,c("time","type","aggregate_consumption","prob_contact_I_weighted","I","R","D","infolag")], id_cols=c("time","type","infolag"), names_from=type, values_from=c(aggregate_consumption, prob_contact_I_weighted, I, R, D), names_sep="_") %>%
                        mutate(total_implied_savings = aggregate_consumption_coordinated - aggregate_consumption_decentralized) %>%
                        mutate(PV_decentralized_losses = (58000/365 - aggregate_consumption_decentralized)*discount_factor^time) %>%
                        mutate(PV_coordinated_losses = (58000/365 -aggregate_consumption_coordinated)*discount_factor^time) %>%
                        # mutate(PV_blanket_1_losses = (58000/365 -aggregate_consumption_blanket_1)*discount_factor^time) %>%
                        mutate(PV_total_implied_savings = (aggregate_consumption_coordinated - aggregate_consumption_decentralized)*discount_factor^time) %>%
                        mutate(total_averted_contacts_with_Is = (prob_contact_I_weighted_decentralized - prob_contact_I_weighted_coordinated) ) %>%
                        mutate(implied_savings_per_averted_I_contact = total_implied_savings/total_averted_contacts_with_Is)

      names(scenarios_wide)

      #big_dfrm = big_dfrm %>% filter(type != "blanket_1")

      losses_and_deaths <- scenarios_wide %>% group_by(infolag) %>%
            summarise(decentralized_losses = round(sum(PV_decentralized_losses),2),
                  decentralized_casesp100k = round(max(I_decentralized + R_decentralized + D_decentralized)*100000,0),
                  coordinated_losses = round(sum(PV_coordinated_losses),2),
                  coordinated_casesp100k = round(max(I_coordinated + R_coordinated + D_coordinated)*100000,0)
                  )

      bar_data <- pivot_longer(losses_and_deaths, cols=c("decentralized_losses", "decentralized_casesp100k", "coordinated_losses", "coordinated_casesp100k"), names_to = c("type","measure"), names_sep="_", values_to = "values")

      print(bar_data, nrow=100)

      bar_summary_nolag <- ggplot((bar_data %>% filter(infolag==unique(infolag)[1])), aes(fill=type, x=reorder(type, -values), y=values)) +
            geom_bar(position="dodge", stat="identity") +
            ggtitle(case1_label) + 
            # facet_wrap(~'measure + infolag', scales = "free_y", labeller = as_labeller(c(casesp100k = "Total cases per 100,000",losses="Individual loss ($/person)")), strip.position = "left") +
            facet_wrap(~measure, labeller = as_labeller(c(casesp100k = "Total cases per 100,000",losses="Individual loss ($/person)")), strip.position = "left") +
           geom_text(aes(label=round(values,3)), position=position_dodge(width=0.9), vjust=-0.25) +
            #scale_fill_manual(values=Mix) +
            ylab(NULL) + xlab("") + labs(fill = "Policy type") +
            theme_classic() +
            theme(strip.background = element_blank(),
                 strip.placement = "outside") +
            guides(fill=FALSE)


      bar_summary_lag1 <- ggplot((bar_data %>% filter(infolag==unique(infolag)[2])), aes(fill=type, x=reorder(type, -values), y=values)) +
            geom_bar(position="dodge", stat="identity") +
            ggtitle(case2_label) + 
            # facet_wrap(~'measure + infolag', scales = "free_y", labeller = as_labeller(c(casesp100k = "Total cases per 100,000",losses="Individual loss ($/person)")), strip.position = "left") +
            facet_wrap(~measure, labeller = as_labeller(c(casesp100k = "Total cases per 100,000",losses="Individual loss ($/person)")), strip.position = "left") +
           geom_text(aes(label=round(values,3)), position=position_dodge(width=0.9), vjust=-0.25) +
            #scale_fill_manual(values=Mix) +
            ylab(NULL) + xlab("") + labs(fill = "Policy type") +
            theme_classic() +
            theme(strip.background = element_blank(),
                 strip.placement = "outside") +
            guides(fill=FALSE)


      bar_summary_lag2 <- ggplot((bar_data %>% filter(infolag==unique(infolag)[3])), aes(fill=type, x=reorder(type, -values), y=values)) +
            geom_bar(position="dodge", stat="identity") +
            ggtitle(case3_label) + 
            # facet_wrap(~'measure + infolag', scales = "free_y", labeller = as_labeller(c(casesp100k = "Total cases per 100,000",losses="Individual loss ($/person)")), strip.position = "left") +
            facet_wrap(~measure, labeller = as_labeller(c(casesp100k = "Total cases per 100,000",losses="Individual loss ($/person)")), strip.position = "left") +
           geom_text(aes(label=round(values,3)), position=position_dodge(width=0.9), vjust=-0.25) +
            #scale_fill_manual(values=Mix) +
            ylab(NULL) + xlab("") + labs(fill = "Policy type") +
            theme_classic() +
            theme(strip.background = element_blank(),
                 strip.placement = "outside")

      recovery_threshold = 1.5 # The "recovery threshold" is the maximum deviation from the initial steady state allowed before declaring "the economy has recovered". A threshold of 1 means "the economy has recovered when the output gap relative to the initial steady state is less than 1%"

      dollar_format(prefix="$",suffix="", big.mark=",", largest_with_cents = 100)

      scenario_summaries = big_dfrm %>% 
                        group_by(infolag,type) %>%
                        summarise(time_to_suppression = time[min(which(Reff<1))] ,
                              time_to_peak = time[which.max(I)] ,
                              infection_peak = max(I)*100000,
                              total_deaths = max(D)*total_population,
                              recession_peak = -min(aggregate_consumption_deviation),
                              time_to_recovery = time[(recovery_threshold + aggregate_consumption_deviation > 0)][which(time>time_to_peak)][1])#,
                              # total_losses = total_population*sum( (58000/365 - aggregate_consumption)*discount_factor^time) )
      scenario_summaries[-c(1,2)] = round(scenario_summaries[,-c(1,2)],2)
      summary_stats = as.data.frame(t(scenario_summaries)) %>% row_to_names(row_number=2, remove_rows_above=FALSE)
      summary_stats

      indx <- sapply(summary_stats, is.factor)
      summary_stats[indx] <- lapply(summary_stats[indx], function(x) as.numeric(as.character(x)))

      summary_stats["infolag",] = paste0(summary_stats["infolag",])
      summary_stats["time_to_suppression",] = paste0(summary_stats["time_to_suppression",]," days")
      summary_stats["time_to_peak",] = paste0(summary_stats["time_to_peak",]," days")
      summary_stats["time_to_recovery",] = paste0(summary_stats["time_to_recovery",]," days")
      summary_stats["infection_peak",] = paste0(summary_stats["infection_peak",]," cases per 100,000")
      summary_stats["recession_peak",] = paste0(summary_stats["recession_peak",],"% contraction")
      # summary_stats["total_losses",] = paste0(dollar(round(as.numeric(summary_stats["total_losses",]),0)))
      summary_stats

      ##### Generate table

      #summary_table <- ggtexttable(summary_stats, rows = c("Time to peak", "Time to suppression", "Peak infection", "Recession trough", "Time to recovery", "Average per capita loss"), theme = ttheme(colnames.style=colnames_style(size=18)))
      summary_table <- ggtexttable(summary_stats, rows = c(label,"Time to suppression", "Time to peak", "Peak infection", "Total deaths", "Recession trough", "Time to recovery"), theme = ttheme(colnames.style=colnames_style(size=18)))
      summary_table

      lag_summaries <- (bar_summary_nolag | bar_summary_lag1 | bar_summary_lag2) / summary_table
      return(lag_summaries)

}

# Function to generate just the recession metrics
generate_just_metrics <- function(dfrm, type="eqm_plan") {

    big_dfrm <- dfrm %>% 
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
                    mutate(daily_new_unavoidable_cases = parms$tau*parms$rho_o*S*I)

                    # print(head(big_dfrm))

    if(length(unique(big_dfrm$type))==3&type=="eqm_plan") {
        scenarios_wide <- pivot_wider(big_dfrm[,c("time","type","aggregate_consumption","prob_contact_I_weighted","I","R","D","infolag")], id_cols=c("time","type","infolag"), names_from=type, values_from=c(aggregate_consumption, prob_contact_I_weighted, I, R, D), names_sep="_") %>%
                        mutate(total_implied_savings = aggregate_consumption_coordinated - aggregate_consumption_decentralized) %>%
                        mutate(PV_decentralized_losses = (58000/365 - aggregate_consumption_decentralized)*discount_factor^time) %>%
                        mutate(PV_lockdown_losses = (58000/365 - aggregate_consumption_lockdown)*discount_factor^time) %>%
                        mutate(PV_coordinated_losses = (58000/365 -aggregate_consumption_coordinated)*discount_factor^time) %>%
                        # mutate(PV_blanket_1_losses = (58000/365 -aggregate_consumption_blanket_1)*discount_factor^time) %>%
                        mutate(PV_total_implied_savings = (aggregate_consumption_coordinated - aggregate_consumption_decentralized)*discount_factor^time) %>%
                        mutate(total_averted_contacts_with_Is = (prob_contact_I_weighted_decentralized - prob_contact_I_weighted_coordinated) ) %>%
                        mutate(implied_savings_per_averted_I_contact = total_implied_savings/total_averted_contacts_with_Is)

        losses_and_deaths <- scenarios_wide %>% group_by(infolag) %>%
            summarise(decentralized_losses = round(sum(PV_decentralized_losses, na.rm=TRUE),0),
                  decentralized_casesp100k = round(max(I_decentralized + R_decentralized + D_decentralized)*100000,0),
                  lockdown_losses = round(sum(PV_lockdown_losses, na.rm=TRUE),0),
                  lockdown_casesp100k = round(max(I_lockdown + R_lockdown + D_lockdown)*100000,0),
                  coordinated_losses = round(sum(PV_coordinated_losses, na.rm=TRUE),0),
                  coordinated_casesp100k = round(max(I_coordinated + R_coordinated + D_coordinated)*100000,0)
                  )

            # print(losses_and_deaths)

        bar_data <- pivot_longer(losses_and_deaths, cols=c("decentralized_losses", "decentralized_casesp100k", "lockdown_losses", "lockdown_casesp100k", "coordinated_losses", "coordinated_casesp100k"), names_to = c("type","measure"), names_sep="_", values_to = "values")
    }

    if(length(unique(big_dfrm$type))==5&type=="eqm_plan") {
        scenarios_wide <- pivot_wider(big_dfrm[,c("time","type","aggregate_consumption","prob_contact_I_weighted","I","R","D","infolag")], id_cols=c("time","type","infolag"), names_from=type, values_from=c(aggregate_consumption, prob_contact_I_weighted, I, R, D), names_sep="_") %>%
                        mutate(PV_decentralized_losses = (58000/365 - aggregate_consumption_decentralized)*discount_factor^time) %>%
                        mutate(PV_gesir_losses = (58000/365 - aggregate_consumption_nocontrol)*discount_factor^time) %>%
                        mutate(PV_lockdown_coupl_losses = (58000/365 - aggregate_consumption_lockdown_beh)*discount_factor^time) %>%
                        mutate(PV_lockdown_gesir_losses = (58000/365 - aggregate_consumption_lockdown_mech)*discount_factor^time) %>%
                        mutate(PV_coordinated_losses = (58000/365 -aggregate_consumption_coordinated)*discount_factor^time)

        losses_and_deaths <- scenarios_wide %>% group_by(infolag) %>%
            summarise(decentralized_losses = round(sum(PV_decentralized_losses, na.rm=TRUE),0),
                  decentralized_casesp100k = round(max(I_decentralized + R_decentralized + D_decentralized)*100000,0),
                  gesir_losses = round(sum(PV_gesir_losses, na.rm=TRUE),0),
                  gesir_casesp100k = round(max(I_nocontrol + R_nocontrol + D_nocontrol)*100000,0),
                  ld.decen_losses = round(sum(PV_lockdown_coupl_losses, na.rm=TRUE),0),
                  ld.decen_casesp100k = round(max(I_lockdown_beh + R_lockdown_beh + D_lockdown_beh)*100000,0),
                  lockdown.gesir_losses = round(sum(PV_lockdown_gesir_losses, na.rm=TRUE),0),
                  lockdown.gesir_casesp100k = round(max(I_lockdown_mech + R_lockdown_mech + D_lockdown_mech)*100000,0),
                  coordinated_losses = round(sum(PV_coordinated_losses, na.rm=TRUE),0),
                  coordinated_casesp100k = round(max(I_coordinated + R_coordinated + D_coordinated)*100000,0)
                  )

            # print(losses_and_deaths)

        bar_data <- pivot_longer(losses_and_deaths, cols=c("decentralized_losses", "decentralized_casesp100k", "gesir_losses", "gesir_casesp100k", "lockdown.decen_losses", "lockdown.decen_casesp100k", "lockdown.gesir_losses", "lockdown.gesir_casesp100k", "coordinated_losses", "coordinated_casesp100k"), names_to = c("type","measure"), names_sep="_", values_to = "values")
    }


        if(type=="flexible") {
        scenarios_wide <- pivot_wider(big_dfrm[,c("time","type","aggregate_consumption","prob_contact_I_weighted","I","R","D","infolag")], id_cols=c("time","type","infolag"), names_from=type, values_from=c(aggregate_consumption, prob_contact_I_weighted, I, R, D), names_sep="_")

        print(colnames(scenarios_wide))

            for(i in seq_along(unique(big_dfrm[,"type"]))) {
                core_name <- unique(big_dfrm[,"type"])[i]
                new_colname <- paste0("PV_",core_name,"_losses")
                cons_colname <- paste0("aggregate_consumption_",core_name)
                print(cons_colname)
                scenarios_wide <- scenarios_wide %>% mutate(!!new_colname := (58000/365 - !!rlang::sym(cons_colname))*discount_factor^time)
            }


            temp <- list() 

            for(i in seq_along(unique(big_dfrm[,"type"]))) {
                core_name <- unique(big_dfrm[,"type"])[i]
                loss_colname <- paste0(core_name,"_losses")
                case_colname <- paste0(core_name,"_casesp100k")

                loss_colname_input <- paste0("PV_",core_name,"_losses")

                I_name <- paste0("I_",core_name)
                R_name <- paste0("R_",core_name)
                D_name <- paste0("D_",core_name)

                temp_wide <- scenarios_wide %>% group_by(infolag) %>% 
                    summarise(!!loss_colname := round(sum( !!rlang::sym(loss_colname_input), na.rm=TRUE )),
                                !!case_colname := round(max(!!rlang::sym(D_name) + !!rlang::sym(I_name) + !!rlang::sym(R_name))*100000 ,0 ) )

                temp[[i]] <- pivot_longer(temp_wide, cols=c(
                    !!rlang::sym(loss_colname), !!rlang::sym(case_colname), 
                    ), names_to = c("type","measure"), names_sep="_", values_to = "values")
            }

            losses_and_deaths <- rbindlist(temp)

            bar_data <- losses_and_deaths
        }

        if(type=="eqm_plan2") {
        scenarios_wide <- pivot_wider(big_dfrm[,c("time","type","aggregate_consumption","prob_contact_I_weighted","I","R","D","infolag")], id_cols=c("time","type","infolag"), names_from=type, values_from=c(aggregate_consumption, prob_contact_I_weighted, I, R, D), names_sep="_") %>%
                        mutate(PV_no.control_losses = (58000/365 - aggregate_consumption_no.control)*discount_factor^time) %>%
                        mutate(PV_voluntary.isolation_losses = (58000/365 - aggregate_consumption_voluntary.isolation)*discount_factor^time) %>%
                        mutate(PV_lockdown_losses = (58000/365 - aggregate_consumption_lockdown)*discount_factor^time) %>%
                        mutate(PV_casemin_losses = (58000/365 - aggregate_consumption_ld.casemin)*discount_factor^time) %>%
                        mutate(PV_targeted.no.casecost_losses = (58000/365 - aggregate_consumption_targeted.no.casecost)*discount_factor^time) %>%
                        mutate(PV_targeted.with.casecost_losses = (58000/365 - aggregate_consumption_targeted.with.casecost)*discount_factor^time)

        losses_and_deaths <- scenarios_wide %>% group_by(infolag) %>%
            summarise(
                no.control_losses = round(sum(PV_no.control_losses, na.rm=TRUE),0),
                  no.control_casesp100k = round(max(I_no.control + R_no.control + D_no.control)*100000,0),
                  voluntary.isolation_losses = round(sum(PV_voluntary.isolation_losses, na.rm=TRUE),0),
                  voluntary.isolation_casesp100k = round(max(I_voluntary.isolation + R_voluntary.isolation + D_voluntary.isolation)*100000,0),
                  lockdown_losses = round(sum(PV_lockdown_losses, na.rm=TRUE),0),
                  lockdown_casesp100k = round(max(I_lockdown + R_lockdown + D_lockdown)*100000,0),
                  casemin_losses = round(sum(PV_casemin_losses, na.rm=TRUE),0),
                  casemin_casesp100k = round(max(I_ld.casemin + R_ld.casemin + D_ld.casemin)*100000,0),
                  targeted.no.casecost_losses = round(sum(PV_targeted.no.casecost_losses, na.rm=TRUE),0),
                  targeted.no.casecost_casesp100k = round(max(I_targeted.no.casecost + R_targeted.no.casecost + D_targeted.no.casecost)*100000,0),
                  targeted.with.casecost_losses = round(sum(PV_targeted.with.casecost_losses, na.rm=TRUE),0),
                  targeted.with.casecost_casesp100k = round(max(I_targeted.with.casecost + R_targeted.with.casecost + D_targeted.with.casecost)*100000,0)
                  )

            # print(losses_and_deaths)

        bar_data <- pivot_longer(losses_and_deaths, cols=c(
            "no.control_losses", "no.control_casesp100k", 
            "voluntary.isolation_losses", "voluntary.isolation_casesp100k", 
            "lockdown_losses", "lockdown_casesp100k", 
            "casemin_losses", "casemin_casesp100k", 
            "targeted.no.casecost_losses", "targeted.no.casecost_casesp100k", 
            "targeted.with.casecost_losses", "targeted.with.casecost_casesp100k"
            ), names_to = c("type","measure"), names_sep="_", values_to = "values")
    }

    return(bar_data)

}

# Function to calculate ratio of decentralized to centralized
calculate_ratio <- function(metrics, metric_type) {
    eqm <- metrics %>% filter(type=="decentralized" & measure==metric_type) %>% select(values) %>% as.numeric()

    plan <- metrics %>% filter(type=="coordinated" & measure==metric_type) %>% select(values) %>% as.numeric()

    ratio <- 1 - plan/eqm

    return(ratio)
}

# Function to compare ratios
compare_ratios <- function(list_of_metrics, metric_type) {
    metrics_vec <- rep(NA,length.out=length(list_of_metrics))
    for(i in seq_along(list_of_metrics)) {
        metrics_vec[i] <- calculate_ratio(list_of_metrics[[i]], metric_type)
    }

    return(metrics_vec)
}


# Function to generate recession metrics figures for information frictions scenarios
generate_recession_metrics_small_composite <- function(dfrm, label, plot_type="eqm_plan", end_time=250) {

      label <- as.character(label)

      dfrm <- dfrm %>% filter(time <= end_time)

    if(plot_type=="case_penalty"){
        dfrm$type <- recode(dfrm$type, 
                            no.control = "\nNo control\n", 
                            penalty.0 = "\nTargeted\nisolation\n",
                            penalty.1000 = "\nTargeted\nisolation\nw/ case penalty")
        cols <- c("\nTargeted\nisolation\n" = viridis(3)[3], "\nTargeted\nisolation\nw/ case penalty" = viridis(3)[1],  "\nNo control\n" = viridis(3)[2])
    }

      bar_data <- generate_just_metrics(dfrm, type="flexible")

      write.csv(bar_data, file=paste0(label,"_summary_stats.csv"))

    if(length(unique(bar_data$type))==3&plot_type!="case_penalty"){
        bar_data$type <- recode(bar_data$type, decentralized = "\nVoluntary\nisolation\n", coordinated = "\nTargeted\nisolation\n", nocontrol = "\nNo control\n")
        cols <- c("\nTargeted\nisolation\n" = viridis(3)[3], "\nVoluntary\nisolation\n" = viridis(3)[1], "\nNo control\n" = viridis(3)[2])
    }

      bar_summary <- ggplot((bar_data %>% filter(infolag==unique(infolag)[1])), aes(fill=type, x=reorder(type, -values), y=values)) +
            geom_bar(position="dodge", stat="identity") +
            facet_wrap(~measure, labeller = as_labeller(c(casesp100k = "Total cases per 100,000",losses="Individual loss ($/person)")), strip.position = "left") +
           geom_text(aes(label=round(values,3)), position=position_dodge(width=0.9), vjust=-0.25) +
            ylab(NULL) + xlab("") + labs(fill = "Policy type") +
            theme_classic() +
            theme(strip.background = element_blank(),
                 strip.placement = "outside",
                 text = element_text(size=24)) + 
            scale_color_viridis() +
            guides(fill=FALSE, color=FALSE)


    if(length(unique(dfrm$type))==3){ bar_summary <- bar_summary + scale_color_manual(values=cols) }

        biggest_loss <- max((bar_data %>% filter(infolag==unique(infolag)[1] & measure=="losses"))$values)

      bar_summary_losses <- ggplot((bar_data %>% filter(infolag==unique(infolag)[1] & measure=="losses")), aes(fill=type, x=reorder(type, -values), y=values)) +
            geom_bar(position="dodge", stat="identity") +
           geom_text(aes(label=round(values,3)), position=position_dodge(width=0.9), vjust=-0.25) +
            ylab(NULL) + xlab("") + labs(fill = "Control type", title = "Average loss per person ($)") +
            theme_classic() +
            ylim(c(0,(biggest_loss + 500))) +
            theme(strip.background = element_blank(),
                 strip.placement = "outside",
                 text = element_text(size=24)) +
            scale_color_viridis() +
            guides(fill=FALSE, color=FALSE)


    if(length(unique(dfrm$type))==3){ bar_summary_losses <- bar_summary_losses + scale_fill_manual(values=cols)  }

        smallest_cases <- min((bar_data %>% filter(infolag==unique(infolag)[1] & measure=="casesp100k"))$values)
        biggest_cases <- max((bar_data %>% filter(infolag==unique(infolag)[1] & measure=="casesp100k"))$values)

      bar_summary_cases <- ggplot((bar_data %>% filter(infolag==unique(infolag)[1] & measure=="casesp100k")), aes(fill=type, x=reorder(type, -values), y=values)) +
            # geom_bar(position="dodge", stat="identity") +
            geom_point(aes(color=type), size=10) +
           geom_text(aes(label=round(values,3)), position=position_dodge(width=0.9), vjust=-1.5) +
            ylab(NULL) + xlab("") + labs(fill = "Control type", color = "Control type", title = "Total cases/100,000") +
            theme_classic() +
            ylim(c( (smallest_cases - 500) , (biggest_cases + 10000) )) +
            theme(strip.background = element_blank(),
                 strip.placement = "outside",
                 text = element_text(size=24)) +
            guides(fill=FALSE, color=FALSE)

    if(length(unique(dfrm$type))==3){ bar_summary_cases <- bar_summary_cases + scale_color_manual(values=cols)  }

      lag_summaries <- list(overall=bar_summary, losses=bar_summary_losses, cases=bar_summary_cases)
      return(lag_summaries)
}

##### CODE TO READ IN ALL CSV FILES IN THE SPECIFIED WORKING DIRECTORY
read_vpfn_csvs <- function(plot_dir) {
  directory <- paste0("../../Results/",plot_dir)
  setwd(directory)

  vpf_filenames = list.files(pattern="*.csv")

  value_policy_function_list = list()
  for(name in seq_along(vpf_filenames)) {
    
    #creates label to be put into stacked set of value functions
    label = gsub(".csv", "", vpf_filenames[name])
    label = gsub("value_policy_functions__", "", label)
    
    #creates set of variables to determine what each value function corresponds to -- will later be used to extract parts to create eventual value function
    vpfn = grepl("vpfn",label,fixed=TRUE)
    
    eqm = grepl("eqm",label,fixed=TRUE)
    planner = grepl("planner",label,fixed=TRUE)
    # First pass: grab all the descriptors in full
    elasticity <- str_match(label, "elasticity_[-\\d]+(?:\\.\\d+)?")
    gdppc <- str_match(label, "gdppc_[-\\d]+(?:\\.\\d+)?")
    laborsupply <- str_match(label, "labor_[-\\d]+(?:\\.\\d+)?")
    cfr <- str_match(label, "cfr_[-\\d]+(?:\\.\\d+)?")
    R0 <- str_match(label, "R0_[-\\d]+(?:\\.\\d+)?")
    # Second pass: strip down to just numbers
    elasticity <- gsub("elasticity_", "", elasticity)
    gdppc <- gsub("gdppc_", "", gdppc)
    laborsupply <- gsub("labor_", "", laborsupply)
    cfr <- gsub("cfr_", "", cfr)
    R0 <- gsub("R0_", "", R0)

    value_policy_function_list[[name]] = cbind(read_csv(paste0(vpf_filenames[name])),
                                               run=paste0(label), elasticity=elasticity, gdppc=gdppc, laborsupply=laborsupply, cfr=cfr, R0=R0, eqm=eqm, planner=planner)
  }

  #Stacks all policy functions together into one large list
  # stacked_solved_values <- rbindlist(value_policy_function_list)

  # return(list(stacked_solved_values,label))
  return(value_policy_function_list)
}
