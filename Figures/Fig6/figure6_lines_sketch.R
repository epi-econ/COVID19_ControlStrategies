######
# Figure 4
######

library(ggpubr)
library(tidyverse)
library(patchwork)
library(data.table)
library(janitor)
library(scales)

source("../../Code/methods_paper/functions.R")

Mix<-c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99")

### Benchmark levels
Cbm <- 58000/365
Lbm <- 0.3333*24
daily_cons_contacts = 5.166426
daily_work_contacts = 7.513051  
daily_other_contacts = 3.548544     

parms <- read.csv("../../Code/methods_paper/calibration_benchmarks.csv")
parms$kappa_c <- 1
parms$kappa_l <- 1

discount_rate = 0.04                                              # annual discount rate
discount_factor = (1/(1+discount_rate))^(1/365)       # daily discount factor
total_population = 331002651

p_M <- 0.199
t_a <- 0.2737

#####
# Make heatmap
#####

setwd("../../Results/Figure3_Sensitivity_analysis/data/ngm_11_rhoc")

ts_filenames = list.files(pattern="*.csv")

time_series_list = list()
for(name in seq_along(ts_filenames)) {
      label = gsub(".csv", "", ts_filenames[name])
      label = gsub("_0_0.+", "", label)
      label = gsub("ngm_", "", label)
      rhoc = str_extract(label, "(?<=rhoc_)\\d+")
      phi = str_extract(label, "(?<=phi_)\\d+[:punct:]\\d+")
      #print(R0)
      #print(Cbm)
      label = gsub("_rhoc_+.+", "", label)
      time_series_list[[name]] = cbind(read_csv(paste0(ts_filenames[name])),
                              type=paste0(label), rhoc=as.numeric(rhoc), phi=as.numeric(phi))
}

ts_dfrm = rbindlist(time_series_list)

## Small diagnostic
# ts_small <- ts_dfrm %>% filter(rhoc==9 & phi==0.1 & type=="planner")
# base <- ggplot(ts_small, aes(x=time, group=type, linetype=type))
# base + geom_point(aes(y=I)) + theme_bw()

# ts_summary_old = ts_dfrm %>% group_by(rhoc,phi) %>%
#                   summarise(total_cases_ratio = sum(weekly_new_cases[type=="eqm"])/sum(weekly_new_cases[type=="planner"]) ,
#                         excess_recession_severity = sum(aggregate_consumption_deviation[type=="eqm"])/sum(aggregate_consumption_deviation[type=="planner"]) ,
#                         excess_recession_depth = min(aggregate_consumption_deviation[type=="eqm"])/min(aggregate_consumption_deviation[type=="planner"]) ,
#                         excess_eqm_laborS_reduction = min(labor_S[type=="eqm"]) - min(labor_S[type=="planner"]) ,
#                         excess_eqm_laborI_supply = min(labor_I[type=="eqm"]) - min(labor_I[type=="planner"]) ,
#                         fO = min(labor_I[type=="eqm"])/min(labor_S[type=="planner"]),
#                         fI = min(labor_I[type=="planner"])/min(labor_S[type=="eqm"]),
#                         fO_over_fI = fO/fI) %>% 
#                   mutate(productivity_loss = (1 - phi)*100) %>%
#                   ungroup(rhoc,phi)

ts_summary = ts_dfrm %>% group_by(rhoc,phi) %>%
            summarise() %>% 
            mutate(productivity_loss = (1 - phi)*100) %>%
            ungroup(rhoc,phi)

ts_summary
ts_summary = ts_summary %>% filter(phi != "NA")
colnames(ts_summary)[1] = "consumption_contacts"
# View(ts_summary)
# ts_summary = ts_summary %>% filter(excess_eqm_laborS_reduction != -Inf)
ts_summary = ts_summary %>% filter(phi != 0.95)

# write_csv(ts_summary, path="../../../Figures_data_FINAL/figure_3_data_OLD.csv")

scenarios_summary = ts_dfrm %>% 
                  group_by(type,rhoc,phi) %>%
                  summarise(PV_losses = sum(((58000/365) - aggregate_consumption)*discount_factor^time), 
                            cases_per_100k = max(I + R + D)*100000)

scenarios_summary_wide =  pivot_wider(scenarios_summary, id_cols=c("type","rhoc","phi"), names_from=type, values_from=c(PV_losses, cases_per_100k), names_sep="_") %>%
            mutate(economic_loss_ratio = 1 - PV_losses_planner/PV_losses_eqm) %>%
            mutate(total_economic_savings = (PV_losses_eqm - PV_losses_planner)*total_population) %>%
            mutate(cases_per_100k_ratio = cases_per_100k_eqm/cases_per_100k_planner) %>%
            mutate(cases_per_100k_difference = cases_per_100k_eqm - cases_per_100k_planner) %>%
            mutate(productivity_loss = (1 - phi)*100) %>%
            mutate(asymptomatic_prop = ((phi - t_a)/(1-t_a))*100 ) %>%
            mutate(c_l_contact_ratio = rhoc/7.513051) %>%
            mutate(avoid_unavoid_contact_ratio = (rhoc+7.513051)/3.548544 )

PV_losses_adjust = scenarios_summary_wide %>% filter(rhoc==5 & phi==0.9) %>% summarise(factor = 0.91/economic_loss_ratio) %>% pull(factor)

PV_losses_adjust

cl_scenarios_summary_wide_adjusted = scenarios_summary_wide %>% mutate(economic_loss_ratio_ADJUSTED = pmin(economic_loss_ratio*PV_losses_adjust,1) )

write_csv(cl_scenarios_summary_wide_adjusted, path="../../../Figures_data_FINAL/figure_4_cl_heatmap_data.csv")
summary(cl_scenarios_summary_wide_adjusted)

# heatmap_base_raw = ggplot(data = cl_scenarios_summary_wide_adjusted, aes(x=productivity_loss, y=rhoc))

#####
# unavoidable-avoidable ratio
#####

setwd("../ngm_10_rhoo")

ts_filenames = list.files(pattern="*.csv")

time_series_list = list()
for(name in seq_along(ts_filenames)) {
      label = gsub(".csv", "", ts_filenames[name])
      label = gsub("_0_0.+", "", label)
      label = gsub("ngm_", "", label)
      rhoo = str_extract(label, "(?<=rhoo_)\\d+")
      phi = str_extract(label, "(?<=phi_)\\d+[:punct:]\\d+")
      #print(R0)
      #print(Cbm)
      label = gsub("_rhoo_+.+", "", label)
      time_series_list[[name]] = cbind(read_csv(paste0(ts_filenames[name])),
                              type=paste0(label), rhoo=as.numeric(rhoo), phi=as.numeric(phi))
}

ts_dfrm = rbindlist(time_series_list)

ts_summary = ts_dfrm %>% group_by(rhoo,phi) %>%
                  summarise(total_cases_ratio = sum(weekly_new_cases[type=="eqm"])/sum(weekly_new_cases[type=="planner"]) ,
                  		eqm_deviation = sum(aggregate_consumption_deviation[type=="eqm"]),
                  		plan_deviation = sum(aggregate_consumption_deviation[type=="planner"]),
                        excess_recession_severity = sum(aggregate_consumption_deviation[type=="eqm"])/sum(aggregate_consumption_deviation[type=="planner"]) ,
                        excess_recession_depth = min(aggregate_consumption_deviation[type=="eqm"])/min(aggregate_consumption_deviation[type=="planner"]) ,
                        excess_eqm_laborS_reduction = min(labor_S[type=="eqm"]) - min(labor_S[type=="planner"]) ,
                        excess_eqm_laborI_supply = min(labor_I[type=="eqm"]) - min(labor_I[type=="planner"]) ,
                        fO = min(labor_I[type=="eqm"])/min(labor_S[type=="planner"]),
                        fI = min(labor_I[type=="planner"])/min(labor_S[type=="eqm"]),
                        fO_over_fI = fO/fI) %>% 
                  mutate(productivity_loss = (1 - phi)*100) %>%
                  ungroup(rhoo,phi)

# ts_summary = ts_dfrm %>% group_by(rhoo,phi) %>%
#             summarise() %>% 
#             mutate(productivity_loss = (1 - phi)*100) %>%
#             ungroup(rhoo,phi)

ts_summary
ts_summary = ts_summary %>% filter(phi != "NA")
# colnames(ts_summary)[1] = "unavoidable_contacts"
# View(ts_summary)
ts_summary = ts_summary %>% filter(excess_eqm_laborS_reduction != -Inf)
ts_summary = ts_summary %>% filter(phi != 0.95)

# write_csv(ts_summary, path="../../../Figures_data_FINAL/figure_4_data_ua.csv")

ua_scenarios_summary = ts_dfrm %>% 
                  group_by(type,rhoo,phi) %>%
                  summarise(PV_losses = sum(((58000/365) - aggregate_consumption)*discount_factor^time), 
                            cases_per_100k = max(I + R + D)*100000)

ua_scenarios_summary_wide =  pivot_wider(ua_scenarios_summary, id_cols=c("type","rhoo","phi"), names_from=type, values_from=c(PV_losses, cases_per_100k), names_sep="_") %>%
            mutate(economic_loss_ratio = 1 - PV_losses_planner/PV_losses_eqm) %>%
            mutate(total_economic_savings = (PV_losses_eqm - PV_losses_planner)*total_population) %>%
            mutate(cases_per_100k_ratio = cases_per_100k_eqm/cases_per_100k_planner) %>%
            mutate(cases_per_100k_difference = cases_per_100k_eqm - cases_per_100k_planner) %>%
            mutate(productivity_loss = (1 - phi)*100) %>%
            mutate(asymptomatic_prop = ((phi - t_a)/(1-t_a))*100 ) %>%
            mutate(unavoid_avoid_contact_ratio = rhoo/(5.166426+7.513051) )

ua_scenarios_summary_wide <- ua_scenarios_summary_wide %>% filter(unavoid_avoid_contact_ratio != 0)

PV_losses_adjust = ua_scenarios_summary_wide %>% filter(rhoo==4 & phi==0.8) %>% summarise(factor = 0.91/economic_loss_ratio) %>% pull(factor)

PV_losses_adjust

ua_scenarios_summary_wide_adjusted = ua_scenarios_summary_wide %>% mutate(economic_loss_ratio_ADJUSTED = pmin(economic_loss_ratio*PV_losses_adjust,1) )

write_csv(ua_scenarios_summary_wide_adjusted, path="../../../Figures_data_FINAL/figure_4_ua_heatmap_data.csv")
summary(ua_scenarios_summary_wide_adjusted)

#####
# make the maps
#####

heatmap_base_cl = ggplot(data = cl_scenarios_summary_wide_adjusted, aes(x=productivity_loss, y=c_l_contact_ratio))

heatmap_base_ua = ggplot(data = ua_scenarios_summary_wide_adjusted, aes(x=productivity_loss, y=unavoid_avoid_contact_ratio))

cases_ratio_cl = heatmap_base_cl + 
            geom_tile(aes(fill = cases_per_100k_ratio)) +
            # geom_contour_filled(aes(z = excess_recession_severity)) +
            geom_point(aes(x=(1-0.8555)*100, y=daily_cons_contacts/daily_work_contacts), color="black", fill="white", size=5, pch=21) + 
            theme_classic() +
            # ylim(c(0,9)) +
            labs(y = "Ratio of consumption to labor contacts", x = "Productivity loss from infection (%)", title = "Cases/100k averted", fill = "Cases/100k averted")

recession_ratio_cl = heatmap_base_cl + 
            geom_tile(aes(fill = economic_loss_ratio_ADJUSTED)) +
            # geom_contour_filled(aes(z = excess_recession_severity)) +
            geom_point(aes(x=(1-0.8555)*100, y=daily_cons_contacts/daily_work_contacts), color="black", fill="white", size=5, pch=21) + 
            theme_classic() +
            # ylim(c(0,9)) +
            labs(y = "Ratio of consumption to labor contacts", x = "Productivity loss from infection (%)", title = "Individual loss averted", fill = "Loss averted (%)")

cases_ratio_au = heatmap_base_ua + 
            geom_tile(aes(fill = cases_per_100k_ratio)) +
            # geom_contour_filled(aes(z = excess_recession_severity)) +
            geom_point(aes(x=(1-0.8555)*100, y=daily_other_contacts/(daily_cons_contacts+daily_work_contacts)), color="black", fill="white", size=5, pch=21) + 
            theme_classic() +
            # ylim(c(0,9)) +
            labs(y = "Ratio of unavoidable to avoidable contacts", x = "Productivity loss from infection (%)", title = "Cases/100k averted", fill = "Cases/100k averted")

recession_ratio_au = heatmap_base_ua + 
            geom_tile(aes(fill = economic_loss_ratio_ADJUSTED)) +
            # geom_contour_filled(aes(z = excess_recession_severity)) +
            geom_point(aes(x=(1-0.8555)*100, y=daily_other_contacts/(daily_cons_contacts+daily_work_contacts)), color="black", fill="white", size=5, pch=21) + 
            theme_classic() +
            # ylim(c(0,9)) +
            labs(y = "Ratio of unavoidable to avoidable contacts", x = "Productivity loss from infection (%)", title = "Individual loss averted", fill = "Loss averted (%)")

cl_panel <- recession_ratio_cl | cases_ratio_cl

ua_panel <- recession_ratio_au | cases_ratio_au

cl_panel / ua_panel

#####
# Compare contact orderings
#####

setwd("../../../Figure4_Analysis_of_contacts/contact_ordering/old")

filenames = list.files(pattern="*.csv")

list_of_files = list()
for(name in seq_along(filenames)) {
      label = gsub(".csv", "", filenames[name])
      label = str_extract(label, ".+?(?=_)") # extract characters up to the first _ . From: https://stackoverflow.com/questions/40113963/how-to-extract-everything-until-first-occurrence-of-pattern
      label_sep = str_split(label, "-")
      type = label_sep[[1]][1]
      contact_ordering = label_sep[[1]][2]
      list_of_files[[name]] = cbind(read_csv(paste0(filenames[name])), 
            type = type,
            contact_ordering = contact_ordering )
}

all_orderings = rbindlist(list_of_files)
all_orderings = all_orderings %>% 
	mutate(weighted_labor_S = S*labor_S) %>%		
	mutate(weighted_labor_I = I*labor_I) %>%
	mutate(weighted_labor_R = R*labor_R) %>% 
 	mutate(type = recode(type, eqm = "decentralized", plan = "coordinated")) %>%
	mutate(total_contacts = consumption_contacts + labor_contacts + other_contacts) %>%
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


scenarios_wide = pivot_wider(all_orderings[,c("time","type","contact_ordering","aggregate_consumption","prob_contact_I_weighted","I","R","D")], id_cols=c("time","type","contact_ordering"), names_from=c("type"), values_from=c(aggregate_consumption, prob_contact_I_weighted, I, R, D), names_sep="_") %>%

	mutate(total_implied_savings = aggregate_consumption_coordinated - aggregate_consumption_decentralized) %>%
	mutate(PV_decentralized_losses = (58000/365 - aggregate_consumption_decentralized)*discount_factor^time) %>%
	mutate(PV_coordinated_losses = (58000/365 -aggregate_consumption_coordinated)*discount_factor^time) %>%
	#mutate(PV_blanket_1_losses_concave = (58000/365 -aggregate_consumption_blanket_1_concave)*discount_factor^time) %>%
	mutate(PV_total_implied_savings = (aggregate_consumption_coordinated - aggregate_consumption_decentralized)*discount_factor^time)

names(scenarios_wide)

bar_data <- scenarios_wide %>% 
	group_by(contact_ordering) %>%
	summarise(
		savings_from_plan = round((1 - sum(PV_coordinated_losses)/sum(PV_decentralized_losses) )*100,1),
		cases_per_100k_averted = round((max(I_decentralized + R_decentralized + D_decentralized)/max(I_coordinated + R_coordinated + D_coordinated)),1)	)
bar_data

PV_losses_adjust = bar_data %>% filter(contact_ordering=="linear") %>% summarise(factor = 91/savings_from_plan) %>% pull(factor)

PV_losses_adjust

bar_data_adjusted = bar_data %>% mutate(savings_from_plan = round( (savings_from_plan*PV_losses_adjust)/100,2))
bar_data_adjusted

bar_data_long <- bar_data_adjusted %>%
	pivot_longer(!contact_ordering, values_to="value", names_to="outcome")

write_csv(bar_data_long, path="../../../Figures_data_FINAL/figure_4_barplot_data.csv")

bar_summary = ggplot(bar_data_long, aes(fill=contact_ordering, x=reorder(contact_ordering, -value), y=value)) +
      geom_bar(position="dodge", stat="identity", width = 0.75) +
      ggtitle("Total cases per 100,000 and individual losses averted\n(targeted isolation relative to voluntary isolation)") + 
      facet_wrap(vars(outcome), scales = "free_y", labeller = as_labeller(c(savings_from_plan="Individual loss averted (%)",cases_per_100k_averted = "Total cases per 100,000 (ratio)")), strip.position = "left") +
     geom_text(aes(label=round(value,3)), position=position_dodge(width=0.9), vjust=-0.25) +
      #scale_fill_manual(values=Mix) +
      ylab(NULL) + xlab("") + labs(fill = "Contact ordering") +
      theme_classic() +
      ylim(c(0,1)) +
      theme(strip.background = element_blank(),
           strip.placement = "outside")
bar_summary

#####
# Contact function
#####    

### grid of consumption and labor

consn <- seq(0,Cbm,length.out=20)
labor <- seq(0,Lbm,length.out=20)
choice_grid <- expand.grid(consn,labor)
colnames(choice_grid)[1:2] <- c("c_S", "l_S")

choice_grid_big <- cbind(choice_grid, c_I = 58000/365, l_I = 0.3333*24)

### Linear contact function
linear_contact_surface_data <- data.frame(choice_grid, contacts=contacts(parms,choice_grid_big))

gg_linear_contact_surface <- ggplot(data = linear_contact_surface_data) + 
      geom_tile(aes(x = c_S, y = l_S, fill = contacts)) + 
      geom_contour(aes(x = c_S, y = l_S, z = contacts), color="black") +
      geom_point(aes(x=0,y=0), color="black", fill="white", size=5, pch=21) + 
      geom_point(aes(x=58000/365,y=0), color="black", fill="white", size=5, pch=21) + 
      geom_point(aes(x=0,y=0.3333*24), color="black", fill="white", size=5, pch=21) + 
      labs(y="Labor supply by susceptible (hrs/day)", x="", title="Linear contact function", fill="") +
      theme_classic()
gg_linear_contact_surface

### Concave contact function
parms_concave <- parms
parms_concave$kappa_c <- 0.5
parms_concave$kappa_l <- 0.5
parms_concave$rho_c <- daily_cons_contacts/(Cbm^(2*parms_concave$kappa_c))   # daily_cons_contacts = rho_c*Cbm*Cbm
parms_concave$rho_l <- daily_work_contacts/(Lbm^(2*parms_concave$kappa_l))   # daily_work_contacts = rho_l*Lbm*Lbm

concave_contact_surface_data <- data.frame(choice_grid, contacts=contacts(parms_concave,choice_grid_big))

gg_concave_contact_surface <- ggplot(data = concave_contact_surface_data) + 
      geom_tile(aes(x = c_S, y = l_S, fill = contacts)) + 
      geom_contour(aes(x = c_S, y = l_S, z = contacts), color="black") +
      geom_point(aes(x=0,y=0), color="black", fill="white", size=5, pch=21) + 
      geom_point(aes(x=58000/365,y=0), color="black", fill="white", size=5, pch=21) + 
      geom_point(aes(x=0,y=0.3333*24), color="black", fill="white", size=5, pch=21) + 
      labs(x="Consumption by susceptible ($/day)", y="", title="Concave contact function", fill="") +
      theme_classic()
gg_concave_contact_surface

### Convex contact function
parms_convex <- parms
parms_convex$kappa_c <- 2
parms_convex$kappa_l <- 2
parms_convex$rho_c <- daily_cons_contacts/(Cbm^(2*parms_convex$kappa_c))   # daily_cons_contacts = rho_c*Cbm*Cbm
parms_convex$rho_l <- daily_work_contacts/(Lbm^(2*parms_convex$kappa_l))   # daily_work_contacts = rho_l*Lbm*Lbm

convex_contact_surface_data <- data.frame(choice_grid, contacts=contacts(parms_convex,choice_grid_big))

gg_convex_contact_surface <- ggplot(data = convex_contact_surface_data) + 
      geom_tile(aes(x = c_S, y = l_S, fill = contacts)) + 
      geom_contour(aes(x = c_S, y = l_S, z = contacts), color="black") +
      geom_point(aes(x=0,y=0), color="black", fill="white", size=5, pch=21) + 
      geom_point(aes(x=58000/365,y=0), color="black", fill="white", size=5, pch=21) + 
      geom_point(aes(x=0,y=0.3333*24), color="black", fill="white", size=5, pch=21) + 
      labs(title="Convex contact function", x="", y="", fill="Contacts/day") +
      theme_classic()
gg_convex_contact_surface

contact_surface_plots <- gg_convex_contact_surface / gg_linear_contact_surface / gg_concave_contact_surface
contact_surface_plots

write_csv(convex_contact_surface_data, path="../../../Figures_data_FINAL/figure_4_convex_surface_data.csv")
write_csv(linear_contact_surface_data, path="../../../Figures_data_FINAL/figure_4_linear_surface_data.csv")
write_csv(concave_contact_surface_data, path="../../../Figures_data_FINAL/figure_4_concave_surface_data.csv")

#####
# Asymptomatics
# Note that p_a_grid reflects both asymptomatics/pre-symptomatics and minor sufferers, who we assume are able to work
#####

p_a_grid <- seq(0,1,length.out=20)
implied_phi_I <- p_a_grid + t_a*(1-p_a_grid)

pa_to_phiI <- data.frame(asymptomatic = p_a_grid, productivity_loss = 1 - pmin(implied_phi_I,1))

asymptomatics_and_productivity <- ggplot(data = pa_to_phiI) + 
      geom_line(aes(x=asymptomatic, y=productivity_loss), size=1) +
      labs(x="Proportion of pre/asymptomatic and minor sufferers", y="Implied productivity loss") +
      theme_bw()
asymptomatics_and_productivity

write_csv(pa_to_phiI, path="../../../Figures_data_FINAL/figure_4_asymptomatic_productivity_mapping.csv")


#####
# Make mockup plot
#####

fig4_mockup <- cl_panel / ua_panel / bar_summary + plot_annotation(tag_levels = "a")

fig4_mockup

png(paste0("../../fig4_sketch.png"), width=1000, height=1000)
fig4_mockup
dev.off()

### Mockup plot with contact surfaces and asymptomatic share

fig4_mockup_2 <- ((asymptomatics_and_productivity / gg_convex_contact_surface / gg_linear_contact_surface / gg_concave_contact_surface) | (recession_ratio_cl | cases_ratio_cl) / (recession_ratio_au | cases_ratio_au) / bar_summary) + plot_annotation(tag_levels = "a")

png(paste0("../../fig4_sketch_with_function_shapes.png"), width=1400, height=1000)
fig4_mockup_2
dev.off()







########################### OLD CODE

#####
# Calculate summary statistics, generate summary table
#####

# with 10^3 points, need to increase recovery threshold
recovery_threshold = 1.5 # The "recovery threshold" is the maximum deviation from the initial steady state allowed before declaring "the economy has recovered". A threshold of 1 means "the economy has recovered when the output gap relative to the initial steady state is less than 1%"

#### Calculate statistics with group_by
##### Equilibrium
scenario_summaries = all_orderings %>% 
                  filter(type == "eqm") %>%
                  group_by(contact_ordering) %>%
                  summarise(time_to_suppression = time[min(which(Reff<1))] ,
                        time_to_peak = time[which.max(I)] ,
                        infection_peak = max(I)*100000,
                        recession_peak = -min(aggregate_consumption_deviation),
                        time_to_recovery = time[(recovery_threshold + aggregate_consumption_deviation > 0)][which(time>time_to_peak)][1],
                        total_losses = total_population*sum( (58000/365 - aggregate_consumption)*discount_factor^time) )
scenario_summaries[-1] = round(scenario_summaries[,-1],2)
summary_stats = as.data.frame(t(scenario_summaries)) %>% row_to_names(row_number=1)
summary_stats
summary_stats <- summary_stats %>% relocate(convex, .after=linear)
summary_stats
colnames(summary_stats) = c("Concave","Convex","Linear")

indx <- sapply(summary_stats, is.factor)
summary_stats[indx] <- lapply(summary_stats[indx], function(x) as.numeric(as.character(x)))

summary_stats["time_to_suppression",] = paste0(summary_stats["time_to_suppression",]," days")
summary_stats["time_to_peak",] = paste0(summary_stats["time_to_peak",]," days")
summary_stats["infection_peak",] = paste0(summary_stats["infection_peak",],"% infected")
summary_stats["recession_peak",] = paste0(summary_stats["recession_peak",],"% contraction")
summary_stats["time_to_recovery",] = paste0(summary_stats["time_to_recovery",]," days")
summary_stats["total_losses",] = paste0(summary_stats["total_losses",]," USD")
colnames(summary_stats) = c("Concave","Convex","Linear")
summary_stats

summary_stats_eqm = summary_stats

##### Planner
scenario_summaries = all_orderings %>% 
                  filter(type == "plan") %>%
                  group_by(contact_ordering) %>%
                  summarise(time_to_suppression = time[min(which(Reff<1))] ,
                        time_to_peak = time[which.max(I)] ,
                        infection_peak = max(I)*100000,
                        recession_peak = -min(aggregate_consumption_deviation),
                        time_to_recovery = time[(recovery_threshold + aggregate_consumption_deviation > 0)][which(time>time_to_peak)][1],
                        total_losses = total_population*sum( (58000/365 - aggregate_consumption)*discount_factor^time) )
scenario_summaries[-1] = round(scenario_summaries[,-1],2)
summary_stats = as.data.frame(t(scenario_summaries)) %>% row_to_names(row_number=1)
summary_stats

indx <- sapply(summary_stats, is.factor)
summary_stats[indx] <- lapply(summary_stats[indx], function(x) as.numeric(as.character(x)))

summary_stats["time_to_suppression",] = paste0(summary_stats["time_to_suppression",]," days")
summary_stats["time_to_peak",] = paste0(summary_stats["time_to_peak",]," days")
summary_stats["infection_peak",] = paste0(summary_stats["infection_peak",],"% infected")
summary_stats["recession_peak",] = paste0(summary_stats["recession_peak",],"% contraction")
summary_stats["time_to_recovery",] = paste0(summary_stats["time_to_recovery",]," days")
summary_stats["total_losses",] = paste0(summary_stats["total_losses",]," USD")
colnames(summary_stats) = c("Concave","Convex","Linear")
summary_stats

summary_stats_plan = summary_stats

##### Generate table

summary_table <- ggtexttable(summary_stats[-1,], rows = c("Time to peak", "Time to suppression", "Peak infection", "Recession trough", "Time to recovery", "Long-run output gap"), 
                        theme = ttheme(colnames.style=colnames_style(size=18)))

summary_table

summary_table_eqm <- ggtexttable(summary_stats_eqm[-1,], rows = c("Time to suppression", "Time to peak", "Peak infection", "Recession trough", "Time to recovery", "Economy-wide loss"), 
                        theme = ttheme(colnames.style=colnames_style(size=18)))
summary_table_eqm

summary_table_plan <- ggtexttable(summary_stats_plan[-1,], rows = c("Time to suppression", "Time to peak", "Peak infection", "Recession trough", "Time to recovery", "Economy-wide loss"), 
                        theme = ttheme(colnames.style=colnames_style(size=18)))


#####
# Write out csvs of data for pretty and final plots
#####

write_csv(summary_stats_plan, path="../../Figures_data_FINAL/figure_4_planner_summary_table.csv",)
write_csv(summary_stats_eqm, path="../../Figures_data_FINAL/figure_4_eqm_summary_table.csv")
write_csv(all_orderings, path="../../Figures_data_FINAL/figure_4_lineplot_data.csv")

#####
# Generate contact orderings plots
#####

eqm_orderings = all_orderings %>% filter(type == "eqm")

eqm_base = ggplot(data=eqm_orderings, aes(x=time))

infection = eqm_base + 
            geom_line(aes(y=I, group=contact_ordering, color=contact_ordering), size=1) +
            geom_line(aes(y=D, group=contact_ordering, color=contact_ordering), linetype="dashed", size=1) +
            theme_bw() + ggtitle("Infecteds (solid) and dead (dashed)") +
            xlab("Day") + ylab("Proportion of population") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

Reff = eqm_base + 
            geom_line(aes(y=Reff, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("R_effective") +
            xlab("Day") + ylab("R_effective") + 
            geom_hline(yintercept=1, linetype="dashed", alpha=0.75, color="darkgray") +
            scale_colour_manual(values=Mix) + guides(color=FALSE)

cases = eqm_base + 
            geom_line(aes(y=weekly_new_cases, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("Daily new cases") +
            xlab("Day") + ylab("Cases") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

infect_prob_S = eqm_base + 
            geom_line(aes(y=prob_infection_S, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("Probability of infection") +
            xlab("Day") + ylab("Probability") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

uninfected = eqm_base + 
            geom_line(aes(y=(S+R), group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("Uninfecteds (effective labor force)") +
            xlab("Day") + ylab("Proportion of population") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

recovereds = eqm_base + 
            geom_hline(yintercept = 1 - 1/all_orderings$Reff[1], linetype = "dashed") +
            geom_line(aes(y=R, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("Total recovered") +
            xlab("Day") + ylab("Proportion of population") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

consumption = eqm_base + 
            geom_line(aes(y=aggregate_consumption_deviation, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("GDP deviation") +
            xlab("Day") + ylab("Percentage of initial steady state") + 
            scale_colour_manual(values=Mix) + labs(linetype="Contact\nordering")

hours = eqm_base + 
            geom_line(aes(y=aggregate_hours_deviation, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("Aggregate hours deviation") +
            xlab("Day") + ylab("Percentage of initial steady state") + 
            scale_colour_manual(values=Mix)

weighted_labor_supplies = eqm_base +
            geom_line(aes(y = S*labor_S, group=contact_ordering, linetype=contact_ordering), size=1, color = "black") +
            geom_line(aes(y = I*labor_I, group=contact_ordering, linetype=contact_ordering), size=1, color = "firebrick4") +
            geom_line(aes(y = R*labor_R, group=contact_ordering, linetype=contact_ordering), size=1, color = "dodgerblue4") +
            theme_bw() + ggtitle("Aggregate labor supplies\n(black = S, red = I, blue = R)") +
            xlab("Day") + ylab("pop_size*hours") + theme(legend.key.size=grid::unit(2,"lines")) + labs(linetype="Contact\nordering")

total_contacts = eqm_base + 
            geom_line(aes(y=total_contacts, group=contact_ordering, color=contact_ordering), size=1) +
            geom_hline(yintercept = 12, linetype="dashed", color="darkgray") +
            geom_hline(yintercept = 5, linetype="dashed", color="darkgray") +
            theme_bw() + ggtitle("Total S-I contacts") +
            xlab("Day") + ylab("Daily contacts") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

labor_S = eqm_base + 
            geom_line(aes(y=labor_S, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("labor supply S") +
            xlab("Day") + ylab("Time spent") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

labor_I = eqm_base + 
            geom_line(aes(y=labor_I, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("labor supply I") +
            xlab("Day") + ylab("Time spent") + ylim(c(2,8)) +
            scale_colour_manual(values=Mix) + guides(color=FALSE)

labor_R = eqm_base + 
            geom_line(aes(y=labor_R, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("labor supply R") +
            xlab("Day") + ylab("Time spent") + ylim(0,max(all_orderings$labor_R)) +
            scale_colour_manual(values=Mix) + guides(color=FALSE)


plot_grid(infection, recovereds, Reff, total_contacts, labor_S, labor_I, consumption, weighted_labor_supplies, summary_table, ncol=3)

plot_grid(infection, recovereds, total_contacts, labor_S, labor_I, hours, ncol=3)

png(paste0("fig4_eqm_contact_orderings_sketch.png"), width=1300, height=900)
plot_grid(infection, recovereds, Reff, total_contacts, labor_S, labor_I, consumption, weighted_labor_supplies, summary_table, ncol=3)
dev.off()

################## planner

plan_orderings = all_orderings %>% filter(type == "plan")

plan_base = ggplot(data=plan_orderings, aes(x=time))

infection = plan_base + 
            geom_line(aes(y=I, group=contact_ordering, color=contact_ordering), size=1) +
            geom_line(aes(y=D, group=contact_ordering, color=contact_ordering), linetype="dashed", size=1) +
            theme_bw() + ggtitle("Infecteds (solid) and dead (dashed)") +
            xlab("Day") + ylab("Proportion of population") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

Reff = plan_base + 
            geom_line(aes(y=Reff, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("R_effective") +
            xlab("Day") + ylab("R_effective") + 
            geom_hline(yintercept=1, linetype="dashed", alpha=0.75, color="darkgray") +
            scale_colour_manual(values=Mix) + guides(color=FALSE)

cases = plan_base + 
            geom_line(aes(y=weekly_new_cases, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("Daily new cases") +
            xlab("Day") + ylab("Cases") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

infect_prob_S = plan_base + 
            geom_line(aes(y=prob_infection_S, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("Probability of infection") +
            xlab("Day") + ylab("Probability") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

uninfected = plan_base + 
            geom_line(aes(y=(S+R), group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("Uninfecteds (effective labor force)") +
            xlab("Day") + ylab("Proportion of population") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

recovereds = plan_base + 
            geom_hline(yintercept = 1 - 1/all_orderings$Reff[1], linetype = "dashed") +
            geom_line(aes(y=R, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("Total recovered") +
            xlab("Day") + ylab("Proportion of population") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

consumption = plan_base + 
            geom_line(aes(y=aggregate_consumption_deviation, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("GDP deviation") +
            xlab("Day") + ylab("Percentage of initial steady state") + 
            scale_colour_manual(values=Mix) + labs(linetype="Contact\nordering")

hours = plan_base + 
            geom_line(aes(y=aggregate_hours_deviation, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("Aggregate hours deviation") +
            xlab("Day") + ylab("Percentage of initial steady state") + 
            scale_colour_manual(values=Mix)

weighted_labor_supplies = plan_base +
            geom_line(aes(y = S*labor_S, group=contact_ordering, linetype=contact_ordering), size=1, color = "black") +
            geom_line(aes(y = I*labor_I, group=contact_ordering, linetype=contact_ordering), size=1, color = "firebrick4") +
            geom_line(aes(y = R*labor_R, group=contact_ordering, linetype=contact_ordering), size=1, color = "dodgerblue4") +
            theme_bw() + ggtitle("Aggregate labor supplies\n(black = S, red = I, blue = R)") +
            xlab("Day") + ylab("pop_size*hours") + theme(legend.key.size=grid::unit(2,"lines")) + labs(linetype="Contact\nordering")

total_contacts = plan_base + 
            geom_line(aes(y=total_contacts, group=contact_ordering, color=contact_ordering), size=1) +
            geom_hline(yintercept = 12, linetype="dashed", color="darkgray") +
            geom_hline(yintercept = 5, linetype="dashed", color="darkgray") +
            theme_bw() + ggtitle("Total S-I contacts") +
            xlab("Day") + ylab("Daily contacts") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE) + theme(text=element_text(family="LM Roman 10", size=15))

labor_S = plan_base + 
            geom_line(aes(y=labor_S, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("labor supply S") +
            xlab("Day") + ylab("Time spent") + ylim(c(0,8)) +
            scale_colour_manual(values=Mix) + guides(color=FALSE)

labor_I = plan_base + 
            geom_line(aes(y=labor_I, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("labor supply I") +
            xlab("Day") + ylab("Time spent") + ylim(c(0,8)) +
            scale_colour_manual(values=Mix) + guides(color=FALSE)

labor_R = plan_base + 
            geom_line(aes(y=labor_R, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("labor supply R") +
            xlab("Day") + ylab("Time spent") + ylim(0,max(all_orderings$labor_R)) +
            scale_colour_manual(values=Mix) + guides(color=FALSE)

plot_grid(infection, recovereds, Reff, total_contacts, labor_S, labor_I, consumption, weighted_labor_supplies, summary_table, ncol=3)

plot_grid(infection, recovereds, total_contacts, labor_S, labor_I, hours, ncol=3)

png(paste0("fig4_plan_contact_orderings_sketch.png"), width=1300, height=900)
plot_grid(infection, recovereds, Reff, total_contacts, labor_S, labor_I, consumption, weighted_labor_supplies, summary_table, ncol=3)
dev.off()


######### Final figure

infection_eqm = eqm_base + 
            geom_line(aes(y=I, group=contact_ordering, color=contact_ordering), size=1) +
            geom_line(aes(y=D, group=contact_ordering, color=contact_ordering), linetype="dashed", size=1) +
            theme_bw() + ggtitle("Infecteds (solid) and dead (dashed)") +
            xlab("Day") + ylab("Proportion of population") + ggtitle("Decentralized epidemic progression") +
            scale_colour_manual(values=Mix) + guides(color=FALSE)

consumption_eqm = eqm_base + 
            geom_line(aes(y=aggregate_consumption_deviation, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("GDP deviation") +
            xlab("Day") + ylab("Percentage change from initial level") + ggtitle("Decentralized economy GDP deviation") +
            scale_colour_manual(values=Mix) + labs(linetype="Contact\nordering")

infection_plan = plan_base + 
            geom_line(aes(y=I, group=contact_ordering, color=contact_ordering), size=1) +
            geom_line(aes(y=D, group=contact_ordering, color=contact_ordering), linetype="dashed", size=1) +
            theme_bw() + ggtitle("Infecteds (solid) and dead (dashed)") +
            xlab("Day") + ylab("Proportion of population") + ggtitle("Planner's epidemic progression") +
            scale_colour_manual(values=Mix) + guides(color=FALSE)

consumption_plan = plan_base + 
            geom_line(aes(y=aggregate_consumption_deviation, group=contact_ordering, color=contact_ordering), size=1) +
            theme_bw() + ggtitle("GDP deviation") +
            xlab("Day") + ylab("Percentage change from initial level") + ggtitle("Planner's economy GDP deviation") +
            scale_colour_manual(values=Mix) + labs(linetype="Contact\nordering")

plot_grid(infection_eqm, consumption_eqm, summary_table_eqm, infection_plan, consumption_plan, summary_table_plan, ncol=3)

png(paste0("fig4.png"), width=1300, height=500)
plot_grid(infection_eqm, consumption_eqm, summary_table_eqm, infection_plan, consumption_plan, summary_table_plan, ncol=3)
dev.off()




#####
# Static plots
#####

setwd("../../Results/Figure4_Analysis_of_contacts/ts_full_big_data")

filenames = list.files(pattern="*.csv")

list_of_files = list()
for(name in seq_along(filenames)) {
      label = gsub(".csv", "", filenames[name])
      benchmark_cons_contacts = str_extract(label, "(?<=rhoc_)\\d+")
      list_of_files[[name]] = cbind(read_csv(paste0(filenames[name])), benchmark_cons_contacts = benchmark_cons_contacts )
}

all_scenarios = rbindlist(list_of_files)
all_scenarios = all_scenarios %>% group_by(benchmark_cons_contacts) %>%
                  mutate(GDPpc_loss = aggregate_consumption[1] - aggregate_consumption) %>%
                  mutate(discount_factor_sequence = discount_factor^time) %>%
                  mutate(PV_GDPpc_loss = discount_factor_sequence*GDPpc_loss) %>%
                  ungroup(benchmark_cons_contacts)

scenarios_summary = all_scenarios %>% group_by(benchmark_cons_contacts) %>% 
                        summarise(total_GDPpc_loss = sum(PV_GDPpc_loss),
                                    avg_daily_GDPpc_loss = total_GDPpc_loss/n(),
                                    total_case_count = sum(weekly_new_cases)) %>%
                        ungroup(benchmark_cons_contacts) %>% 
                        sapply(as.numeric) %>% as.data.frame() %>%
                        mutate(dTotalCase_dConsContacts = c(NA,diff(total_case_count))) %>%
                        mutate(marginal_avg_daily_GDPpc_loss = c(NA,diff(avg_daily_GDPpc_loss)))
scenarios_summary

plots_base = ggplot(data=scenarios_summary, aes(x=benchmark_cons_contacts))

TC_of_contacts = plots_base + geom_line(aes(y=avg_daily_GDPpc_loss), size = 1) + 
      geom_vline(xintercept = 4.2877, linetype= "dashed") +
      labs(x="Average daily consumption activity contacts (#)", y="Annuitized daily GDPpc loss ($)", title="GDPpc loss due to virus and\ninfection-avoiding behavior") + 
      theme_classic()

MC_of_contacts = plots_base + geom_line(aes(y=marginal_avg_daily_GDPpc_loss), size = 1) + 
      geom_vline(xintercept = 4.2877, linetype= "dashed") +
      labs(x="Average daily consumption activity contacts (#)", y="Marginal annuitized daily GDPpc loss ($)", title="Marginal cost of contacts at consumption activities\ngiven infection-avoiding behavior") + 
      theme_classic()

total_new_cases_from_contacts = plots_base + geom_line(aes(y=total_case_count), size = 1) +
      geom_vline(xintercept = 4.2877, linetype= "dashed") +
      labs(x="Average daily consumption activity contacts (#)", y="Total cases (% of population)", title="Total cases given infection-avoiding behavior") + 
      theme_classic()

marginal_new_cases_from_contacts = plots_base + geom_line(aes(y=dTotalCase_dConsContacts), size = 1) +
      geom_vline(xintercept = 4.2877, linetype= "dashed") +
      labs(x="Average daily consumption activity contacts (#)", y="Marginal cases (% of population)", title="Additional cases given infection-avoiding behavior") + 
      theme_classic()

contacts_upperrow <- plot_grid(TC_of_contacts, total_new_cases_from_contacts, MC_of_contacts, marginal_new_cases_from_contacts, nrow=2)

#####
# Time series plots
#####

setwd("~/Dropbox/Corona_epicon/Results/Figure4_Analysis_of_contacts/mainrun_ts_full_big_data")

filenames = list.files(pattern="*.csv")

list_of_files = list()
for(name in seq_along(filenames)) {
      label = gsub(".csv", "", filenames[name])
      list_of_files[[name]] = cbind(read_csv(paste0(filenames[name])), type = label )
}

ts_full = rbindlist(list_of_files)

window_size = 30
eqm_rows = nrow(ts_full[type=="eqm",])
planner_rows = nrow(ts_full[type=="planner",])

eqm_CCEO = rep(NA,length=eqm_rows)
eqm_IEO = rep(NA,length=eqm_rows)

for(row in (window_size+1):eqm_rows) {
      eqm_CCEO[row] = coef(lm(log(aggregate_consumption) ~ log(consumption_contacts), data=ts_full[type=="eqm",][c((row-window_size):row)]))[[2]]
      eqm_IEO[row] = coef(lm(log(aggregate_consumption) ~ log(I), data=ts_full[type=="eqm",][c((row-window_size):row)]))[[2]]
}

planner_CCEO = rep(NA,length=planner_rows)
planner_IEO = rep(NA,length=planner_rows)

for(row in (window_size+1):planner_rows) {
      planner_CCEO[row] = coef(lm(log(aggregate_consumption) ~ log(consumption_contacts), data=ts_full[type=="planner",][c((row-window_size):row)]))[[2]]
      planner_IEO[row] = coef(lm(log(aggregate_consumption) ~ log(I), data=ts_full[type=="planner",][c((row-window_size):row)]))[[2]]
}

CCEO = c(eqm_CCEO,planner_CCEO)
IEO = c(eqm_IEO,planner_IEO)

ts_full_big = data.frame(ts_full,CCEO=CCEO,IEO=IEO)
markerstart = 35
markerend = 100
ts_full_big$marker = ifelse(ts_full_big$time>markerstart&ts_full_big$time<markerend, ts_full_big$marker <- 1, ts_full_big$marker <- 0)

ts_full_short = ts_full_big[which(ts_full_big$time<150),]

labels = as.factor(c("Decentralized", "Coordinated"))
colors = c(Mix[1], Mix[2])

CCEO_full = ggplot(data = ts_full_big, aes(x=time, group = as.factor(type), color = as.factor(type))) +
            geom_line(aes(y=CCEO), size = 1) +
            geom_vline(aes(xintercept = markerstart), linetype = "dashed", alpha = 0.5) +
            geom_vline(aes(xintercept = markerend), linetype = "dashed", alpha = 0.5) +
            labs(x = "time", y = "Elasticity", color = "Behavior") +
            scale_color_manual(values = c(Mix[2], Mix[1])) + guides(color=FALSE) +
            theme_bw() + ggtitle(paste0("Consumption contacts elasticity of per-capita GDP\n(",window_size,"-day rolling window)")) +
            theme(text = element_text(size=20)) 

eqm_CCEO_zoomedin = ggplot(data = ts_full_big[intersect(which(ts_full_big$type=="eqm"),which(ts_full_big$marker==1)),], aes(x=time, group = as.factor(type), color = as.factor(type))) +
            geom_line(aes(y=CCEO), size = 1, color = Mix[2]) +
            geom_vline(aes(xintercept = markerstart), linetype = "dashed", alpha = 0.5) +
            geom_vline(aes(xintercept = markerend), linetype = "dashed", alpha = 0.5) +
            labs(x = "time", y = "Elasticity", color = "Behavior") +
            scale_color_manual(values = c(Mix[2], Mix[1])) + guides(color=FALSE) +
            theme_bw() + ggtitle("Decentralized") +
            theme(text = element_text(size=20))   

planner_CCEO_zoomedin = ggplot(data = ts_full_big[intersect(which(ts_full_big$type=="planner"),which(ts_full_big$marker==1)),], aes(x=time, group = as.factor(type), color = as.factor(type))) +
            geom_line(aes(y=CCEO), size = 1, color = Mix[1]) +
            geom_vline(aes(xintercept = markerstart), linetype = "dashed", alpha = 0.5) +
            geom_vline(aes(xintercept = markerend), linetype = "dashed", alpha = 0.5) +
            labs(x = "time", y = "Elasticity", color = "Behavior") +
            scale_color_manual(values = c(Mix[2], Mix[1])) + guides(color=FALSE) +
            theme_bw() + ggtitle("Centralized") +
            theme(text = element_text(size=20))       

infection_full = ggplot(data = ts_full_big, aes(x=time, group = type, color = type)) +
                  geom_line(aes(y=I), size=1) +
                  geom_line(aes(y=D), size=1, linetype="dotted") +
                  geom_vline(aes(xintercept = markerstart), linetype = "dashed", alpha = 0.5) +
                  geom_vline(aes(xintercept = markerend), linetype = "dashed", alpha = 0.5) +
                  labs(x = "Days", y = "Fraction of population", color = "Behavior") +
                  scale_color_manual(values = c(Mix[2], Mix[1]), labels = labels)  +
                  theme_bw() + ggtitle("Disease progression:\ninfected (solid) and dead (dotted)") +
                  theme(text = element_text(size=20))

recession_full = ggplot(data = ts_full_big, aes(x=time, group = type, color = type)) +
                  geom_line(aes(y=aggregate_consumption_deviation), size=1) +
                  geom_vline(aes(xintercept = markerstart), linetype = "dashed", alpha = 0.5) +
                  geom_vline(aes(xintercept = markerend), linetype = "dashed", alpha = 0.5) +
                  labs(x = "Days", y = "% deviation from initial level") + guides(color = FALSE) +
                  scale_color_manual(values = c(Mix[2], Mix[1]), labels = labels) +
                  theme_bw() + ggtitle("Recession progression:\nchange in GDP") +
                  theme(text = element_text(size=20))

IEO_full = ggplot(data = ts_full_big, aes(x=time, group = as.factor(type), color = as.factor(type))) +
            geom_line(aes(y=IEO), size = 1) +
            geom_vline(aes(xintercept = markerstart), linetype = "dashed", alpha = 0.5) +
            geom_vline(aes(xintercept = markerend), linetype = "dashed", alpha = 0.5) +
            labs(x = "time", y = "Elasticity", color = "Behavior") +
            scale_color_manual(values = c(Mix[2], Mix[1])) + guides(color=FALSE) +
            theme_bw() + ggtitle(paste0("Infection elasticity of per-capita GDP\n(",window_size,"-day rolling window)")) +
            theme(text = element_text(size=20)) 

eqm_IEO_zoomedin = ggplot(data = ts_full_big[intersect(which(ts_full_big$type=="eqm"),which(ts_full_big$marker==1)),], aes(x=time)) +
            geom_line(aes(y=IEO), size = 1, color = Mix[2]) +
            geom_vline(aes(xintercept = markerstart), linetype = "dashed", alpha = 0.5) +
            geom_vline(aes(xintercept = markerend), linetype = "dashed", alpha = 0.5) +
            labs(x = "time", y = "Elasticity", color = "Behavior") +
            #scale_color_manual(values = c(Mix[2], Mix[1])) + guides(color = FALSE) +
            theme_bw() + ggtitle("Decentralized") +
            theme(text = element_text(size=20))

planner_IEO_zoomedin = ggplot(data = ts_full_big[intersect(which(ts_full_big$type=="planner"),which(ts_full_big$marker==1)),], aes(x=time)) +
            geom_line(aes(y=IEO), size = 1, color = Mix[1]) +
            geom_vline(aes(xintercept = markerstart), linetype = "dashed", alpha = 0.5) +
            geom_vline(aes(xintercept = markerend), linetype = "dashed", alpha = 0.5) +
            labs(x = "time", y = "Elasticity", color = "Behavior") +
            #scale_color_manual(values = c(Mix[2], Mix[1])) + guides(color = FALSE) +
            theme_bw() + ggtitle("Centralized") +
            theme(text = element_text(size=20))

weighted_labor_supplies <- ggplot(data = ts_full_big, aes(x=time, group = as.factor(type), linetype = as.factor(type))) +
            geom_line(aes(y=labor_S*S), color="black", size=1) +
            geom_line(aes(y=labor_I*I), color="firebrick4", size=1) +
            geom_line(aes(y=labor_R*R), color="dodgerblue3", size=1) +
            geom_vline(aes(xintercept = markerstart), linetype = "dashed", alpha = 0.5) +
            geom_vline(aes(xintercept = markerend), linetype = "dashed", alpha = 0.5) +
            labs(x = "Days", y = "Average hours worked", linetype = "Behavior") +
            scale_linetype_manual(values = c("solid","dotted"),labels = labels) +
            theme_bw() + ggtitle("Population-weighted labor supplies\n(S = black, I = red, R = blue)") +
            theme(text = element_text(size=20)) 

weighted_labor_supplies_zoomedin <- ggplot(data = ts_full_big[which(ts_full_big$marker==1),], aes(x=time, group = as.factor(type), linetype = as.factor(type))) +
            geom_line(aes(y=labor_S*S), color="black", size=1) +
            geom_line(aes(y=labor_I*I), color="firebrick4", size=1) +
            geom_line(aes(y=labor_R*R), color="dodgerblue3", size=1) +
            geom_vline(aes(xintercept = markerstart), linetype = "dashed", alpha = 0.5) +
            geom_vline(aes(xintercept = markerend), linetype = "dashed", alpha = 0.5) +
            labs(x = "Days", y = "Average hours worked", linetype = "Behavior") +
            scale_linetype_manual(values = c("solid","dotted"),labels = labels) +
            theme_bw() + ggtitle("Population-weighted labor supplies\n(S = black, I = red, R = blue)") +
            theme(text = element_text(size=20)) 

aggregate_outcomes <- plot_grid(infection_full, recession_full, ncol=2)
IEO_zoomedin <- plot_grid(eqm_IEO_zoomedin, planner_IEO_zoomedin, ncol=2, align="v")
IEO_breakdown <- plot_grid(IEO_full, IEO_zoomedin, nrow=2, align = "v")
CCEO_zoomedin <- plot_grid(eqm_CCEO_zoomedin, planner_CCEO_zoomedin, ncol=2, align="v")
CCEO_breakdown <- plot_grid(CCEO_full, CCEO_zoomedin, nrow=2, align = "v")

lower_row <- plot_grid(weighted_labor_supplies, IEO_breakdown, nrow = 1, align = "v")

I_CC_e_GDP <- plot_grid(aggregate_outcomes,lower_row, nrow=2, rel_heights = c(0.4,0.6), align = "v")

#####
# Compiled plots
#####

contacts_analysis <- plot_grid(contacts_upperrow, CCEO_breakdown, ncol=2, align="h")

I_CC_e_GDP

png(file="../contacts_analysis_sketch.png", width = 1200, height = 800)
contacts_analysis
dev.off()
