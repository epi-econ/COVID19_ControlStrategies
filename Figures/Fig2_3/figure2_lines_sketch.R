# Script to make some plots for comparing model features, and data summaries for making final figures.

library(ggpubr)
library(tidyverse)
library(cowplot)
library(patchwork)
library(data.table)
library(janitor)
library(scales)

setwd("../../Code/methods_paper/")

source("functions.R")

Mix<-c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99") 

parms <- read.csv("calibration_benchmarks.csv")
parms$kappa_c <- 1
parms$kappa_l <- 1

discount_rate = 0.04                                              # annual discount rate
discount_factor = (1/(1+discount_rate))^(1/365)       # daily discount factor
total_population = 331002651

setwd("../../Results/Figure2_Models_epi/Data_For_plots")

filenames = list.files(pattern="*.csv")
#filenames = filenames[-1]
#filenames = filenames[c(1,6)]
#filenames = filenames[c(1,2)]

list_of_files = list()
for(name in seq_along(filenames)) {
      label = gsub(".csv", "", filenames[name])
      list_of_files[[name]] = cbind(read.csv(paste0(filenames[name])),type=paste0(label))
}

all_scenarios = rbindlist(list_of_files, fill=TRUE)
all_scenarios = all_scenarios %>% 
                  mutate(weighted_labor_S = S*labor_S) %>%
                  mutate(weighted_labor_I = I*labor_I) %>%
                  mutate(weighted_labor_R = R*labor_R) %>% 
                  mutate(total_contacts = consumption_contacts + labor_contacts + other_contacts) %>%
                  mutate(type = recode(type, eqm_35pt = "decentralized", planner_35pt = "coordinated")) %>%
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

scenarios_wide = pivot_wider(all_scenarios[,c("time","type","aggregate_consumption","prob_contact_I_weighted","I","R","D")], id_cols=c("time","type"), names_from=type, values_from=c(aggregate_consumption, prob_contact_I_weighted, I, R, D), names_sep="_") %>%
                  mutate(total_implied_savings = aggregate_consumption_coordinated - aggregate_consumption_decentralized) %>%
                  mutate(PV_decentralized_losses = (58000/365 - aggregate_consumption_decentralized)*discount_factor^time) %>%
                  mutate(PV_coordinated_losses = (58000/365 -aggregate_consumption_coordinated)*discount_factor^time) %>%
                  mutate(PV_blanket_1_losses = (58000/365 -aggregate_consumption_blanket_1)*discount_factor^time) %>%
                  mutate(PV_total_implied_savings = (aggregate_consumption_coordinated - aggregate_consumption_decentralized)*discount_factor^time) %>%
                  mutate(total_averted_contacts_with_Is = (prob_contact_I_weighted_decentralized - prob_contact_I_weighted_coordinated) ) %>%
                  mutate(implied_savings_per_averted_I_contact = total_implied_savings/total_averted_contacts_with_Is)

names(scenarios_wide)

#all_scenarios = all_scenarios %>% filter(type != "blanket_1")

decentralized_losses = round(sum(scenarios_wide$PV_decentralized_losses),2)
decentralized_deaths = round(max(scenarios_wide$I_decentralized + scenarios_wide$R_decentralized + scenarios_wide$D_decentralized)*100000,0)
coordinated_losses = round(sum(scenarios_wide$PV_coordinated_losses),2)
coordinated_deaths = round(max(scenarios_wide$I_coordinated + scenarios_wide$R_coordinated + scenarios_wide$D_coordinated)*100000,0)
blanket1_losses = round(sum(scenarios_wide$PV_blanket_1_losses),2)
blanket1_deaths = round(max(scenarios_wide$I_blanket_1 + scenarios_wide$R_blanket_1 + scenarios_wide$D_blanket_1)*100000,0)

losses_dfrm = data.frame(blanket = blanket1_losses, decentralized = decentralized_losses, coordinated = coordinated_losses)

epi_dfrm = data.frame(blanket = blanket1_deaths, decentralized = decentralized_deaths, coordinated = coordinated_deaths)

type = c(rep("blanket",2),rep("decentralized",2),rep("coordinated",2))
outcome = rep(c("Loss","CumCases_per_100k"),3)
value = c(blanket1_losses, blanket1_deaths, decentralized_losses, decentralized_deaths, coordinated_losses, coordinated_deaths)
bar_data = data.frame(type,outcome,value)

bar_summary = ggplot(bar_data, aes(fill=type, x=reorder(type, -value), y=value)) +
      geom_bar(position="dodge", stat="identity") +
      ggtitle("Total cases per 100,000 and individual losses") + 
      facet_wrap(~outcome, scales = "free_y", labeller = as_labeller(c(CumCases_per_100k = "Total cases per 100,000",Loss="Individual loss ($/person)")), strip.position = "left") +
     geom_text(aes(label=round(value,3)), position=position_dodge(width=0.9), vjust=-0.25) +
      #scale_fill_manual(values=Mix) +
      ylab(NULL) + xlab("") + labs(fill = "Policy type") +
      theme_classic() +
      theme(strip.background = element_blank(),
           strip.placement = "outside")
bar_summary

#all_scenarios = all_scenarios[which(all_scenarios$time<=100),]

# all_scenarios = all_scenarios %>% group_by(type) %>%
#                   mutate(labor_recovery_rate = diff(aggregate_hours_deviation)/aggregate_hours_deviation) %>% ungroup(type) %>%
#                   filter(time > 3)

#ts_growthrates = ts_growthrates[3:350,]
# head(as.data.frame(ts_growthrates))

# summary(ts_growthrates$labor_recovery_rate)
# which.min(ts_growthrates$labor_recovery_rate)
# which.max(ts_growthrates$labor_recovery_rate)

#####
# Calculate summary statistics, generate summary table
#####

recovery_threshold = 1.5 # The "recovery threshold" is the maximum deviation from the initial steady state allowed before declaring "the economy has recovered". A threshold of 1 means "the economy has recovered when the output gap relative to the initial steady state is less than 1%"

dollar_format(prefix="$",suffix="", big.mark=",", largest_with_cents = 100)

#### Calculate statistics
eqm_summary_stats = all_scenarios %>% 
                  filter(type == "decentralized") %>%
                  summarise(type = "decentralized",
                        time_to_suppression = time[min(which(Reff<1))],
                        time_to_peak = time[which.max(all_scenarios$I)] ,
                        infection_peak = max(I)*100000,
                        recession_peak = -min(aggregate_consumption_deviation),
                        time_to_recovery = time[(recovery_threshold + aggregate_consumption_deviation > 0)][which(time>time_to_peak)][1],
                        total_losses =  decentralized_losses*total_population)  
eqm_summary_stats[,-1] = round(eqm_summary_stats[,-1],2)

plan_summary_stats = all_scenarios %>% 
                  filter(type == "coordinated") %>%
                  summarise(type = "coordinated",
                        time_to_suppression = time[min(which(Reff<1))],
                        time_to_peak = time[which.max(all_scenarios$I)] ,
                        infection_peak = max(I)*100000,
                        recession_peak = -min(aggregate_consumption_deviation),
                        time_to_recovery = time[(recovery_threshold + aggregate_consumption_deviation > 0)][which(time>time_to_peak)][1],
                        total_losses = coordinated_losses*total_population )  
plan_summary_stats[,-1] = round(plan_summary_stats[,-1],2)

blanket_summary_stats = all_scenarios %>% 
                  filter(type == "blanket_1") %>%
                  summarise(type = "blanket_1",
                        time_to_suppression = time[min(which(Reff<1))] ,
                        time_to_peak = time[which.max(all_scenarios$I)] ,
                        infection_peak = max(I)*100000,
                        recession_peak = -min(aggregate_consumption_deviation),
                        time_to_recovery = time[(recovery_threshold + aggregate_consumption_deviation > 0)][which(time>time_to_peak)][1],
                        total_losses = blanket1_losses*total_population )  


blanket_summary_stats[,-1] = round(blanket_summary_stats[,-1],2)

summary_stats = as.data.frame(t(rbind(blanket_summary_stats,eqm_summary_stats, plan_summary_stats))) %>%
                  row_to_names(row_number = 1)
summary_stats

scenario_summaries = all_scenarios %>% 
                  group_by(type) %>%
                  summarise(time_to_suppression = time[min(which(Reff<1))] ,
                        time_to_peak = time[which.max(I)] ,
                        infection_peak = max(I)*100000,
                        total_deaths = max(D)*total_population,
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
summary_stats["time_to_recovery",] = paste0(summary_stats["time_to_recovery",]," days")
summary_stats["infection_peak",] = paste0(summary_stats["infection_peak",]," cases per 100,000")
summary_stats["recession_peak",] = paste0(summary_stats["recession_peak",],"% contraction")
summary_stats["total_losses",] = paste0(dollar(round(as.numeric(summary_stats["total_losses",]),0)))
summary_stats

##### Generate table

#summary_table <- ggtexttable(summary_stats, rows = c("Time to peak", "Time to suppression", "Peak infection", "Recession trough", "Time to recovery", "Average per capita loss"), theme = ttheme(colnames.style=colnames_style(size=18)))
summary_table <- ggtexttable(summary_stats, rows = c("Time to suppression", "Time to peak", "Peak infection", "Recession trough", "Time to recovery", "Total loss"), theme = ttheme(colnames.style=colnames_style(size=18)))
summary_table

#####
# Write out data
#####

all_scenarios_small = all_scenarios[,c("I","weekly_new_cases","aggregate_consumption_deviation","total_contacts","weighted_labor_S","weighted_labor_I","weighted_labor_R","prob_contact_I_weighted","type")] %>% rename(daily_new_cases = weekly_new_cases)

write_csv(all_scenarios_small, path="../../Figures_data_FINAL/figure_2_lineplot_data.csv")
write_csv(bar_data, path="../../Figures_data_FINAL/figure_2_barplot_data.csv")
write_csv(summary_stats, path="../../Figures_data_FINAL/figure_2_summarystats_data.csv")

#####
# Main plots
#####

all_base = ggplot(data=(all_scenarios %>% filter(type!="blanket_1") ), aes(x=time))

infection = all_base + 
            geom_line(aes(y=I*100, group=type, linetype=type, size=type)) +
            #geom_line(aes(y=D, group=type, color=type), linetype="dashed", size=1) +
            theme_bw() + ggtitle("Infecteds") +
            xlab("Day") + ylab("Proportion of population (%)") + 
            scale_size_discrete(range = c(1,1.5)) +
            #scale_colour_manual(values=Mix) + 
            guides(linetype=FALSE, size=FALSE)

infection_zoomed = ggplot(data=(all_scenarios %>% filter(type!="blanket_1", time<75, time>45) ), aes(x=time)) +
            geom_line(aes(y=I*100, group=type, linetype=type, size=type)) +
            #geom_line(aes(y=D, group=type, color=type), linetype="dashed", size=1) +
            theme_bw() + ggtitle("Infecteds: Day 45-75") +
            xlab("Day") + ylab("Proportion of population (%)") + 
            scale_size_discrete(range = c(1,1.5)) +
            #scale_colour_manual(values=Mix) + 
            guides(linetype=FALSE, size=FALSE)

png(paste0("infection_day45-75.png"), width=400, height=300)
infection_zoomed
dev.off()

Reff = all_base + 
            geom_line(aes(y=Reff, group=type, color=type), size=1) +
            theme_bw() + ggtitle("R_effective") +
            xlab("Day") + ylab("R_effective") + 
            geom_hline(yintercept=1, linetype="dashed", alpha=0.75, color="darkgray") +
            scale_colour_manual(values=Mix) + guides(color=FALSE)

cases = all_base + 
            geom_line(aes(y=weekly_new_cases, group=type, color=type), size=1) +
            theme_bw() + ggtitle("Daily new cases") +
            xlab("Day") + ylab("Cases") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

infect_prob_S = all_base + 
            geom_line(aes(y=prob_infection_S, group=type, color=type), size=1) +
            theme_bw() + ggtitle("Probability of infection") +
            xlab("Day") + ylab("Probability") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

uninfected = all_base + 
            geom_line(aes(y=(S+R), group=type, color=type), size=1) +
            theme_bw() + ggtitle("Uninfecteds (effective labor force)") +
            xlab("Day") + ylab("Proportion of population") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

recovereds = all_base + 
            geom_hline(yintercept = 1 - 1/all_scenarios$Reff[1], linetype = "dashed") +
            geom_line(aes(y=R, group=type, color=type), size=1) +
            theme_bw() + ggtitle("Total recovered") +
            xlab("Day") + ylab("Proportion of population") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

consumption = all_base + 
            geom_line(aes(y=aggregate_consumption_deviation, group=type, linetype=type), size=1) +
            theme_bw() + ggtitle("GDP deviation") +
            xlab("Day") + ylab("Percentage of initial steady state") #+ 
            #scale_colour_manual(values=Mix)
hours = all_base + 
            geom_line(aes(y=aggregate_hours_deviation, group=type, color=type), size=1) +
            theme_bw() + ggtitle("Aggregate hours deviation") +
            xlab("Day") + ylab("Percentage of initial steady state") + 
            scale_colour_manual(values=Mix)

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

labor_S = all_base + 
            geom_line(aes(y=labor_S, group=type, color=type), size=1) +
            theme_bw() + ggtitle("labor supply S") +
            xlab("Day") + ylab("Time spent") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

labor_I = all_base + 
            geom_line(aes(y=labor_I, group=type, color=type), size=1) +
            theme_bw() + ggtitle("labor supply I") +
            xlab("Day") + ylab("Time spent") + 
            scale_colour_manual(values=Mix) + guides(color=FALSE)

labor_R = all_base + 
            geom_line(aes(y=labor_R, group=type, color=type), size=1) +
            theme_bw() + ggtitle("labor supply R") +
            xlab("Day") + ylab("Time spent") + ylim(0,max(all_scenarios$labor_R)) +
            scale_colour_manual(values=Mix) + guides(color=FALSE)

prob_contact_I = all_base + 
            geom_line(aes(y=prob_contact_I, group=type, linetype=type), size=1) +
            theme_bw() + ggtitle("Probability given S or R individual contacts a given I-type") +
            xlab("Day") + ylab("Probability (not pop-weighted)") + ylim(0,1)
prob_contact_I


prob_contact_I_wt = all_base + 
            geom_line(aes(y=prob_contact_I_weighted, group=type, linetype=type), size=1) +
            theme_bw() + ggtitle("Probability random S or R contacts any I") +
            xlab("Day") + ylab("Probability") +
            scale_colour_manual(values=Mix) + guides(linetype=FALSE)
prob_contact_I_wt

png(paste0("targeting.png"), width=800, height=500)
plot_grid(prob_contact_I, prob_contact_I_wt)
dev.off()

fig3_mockup <- (weighted_labor_supplies | cases_by_site | contacts_by_site ) / (prob_contact_I_wt | total_contacts) + plot_annotation(tag_levels = "i")

png(paste0("fig3_mockup.png"), width=1000, height=500)
fig3_mockup
dev.off()


fig3_mockup2 <- (weighted_labor_supplies | prob_contact_I_wt ) / (total_contacts | contacts_by_site | cases_by_site ) + plot_annotation(tag_levels = "i")

png(paste0("fig3_mockup2.png"), width=1000, height=500)
fig3_mockup2
dev.off()

# plot_grid(infection, Reff, infect_prob_S, recovereds, total_contacts, hours, labor_S, labor_I, labor_R, ncol=3)

plot_grid(infection, recovereds, Reff, total_contacts, labor_S, labor_I, consumption, weighted_labor_supplies, summary_table, ncol=3)

plot_grid(infection, recovereds, total_contacts, labor_S, labor_I, hours, ncol=3)

mixing_subplot = plot_grid(total_contacts, prob_contact_I_wt, nrow=2, labels = c("d","e"))

png(paste0("fig2_sketch_proper.png"), width=1700, height=750)
plot_grid(infection, consumption, summary_table, mixing_subplot, weighted_labor_supplies, bar_summary, ncol=3, rel_widths = c(1/4,3/8,3/8), labels = c("a","b","c","","f","g"))
dev.off()

png(paste0("fig2_sketch.png"), width=1200, height=900)
plot_grid(infection, recovereds, Reff, total_contacts, labor_S, labor_I, consumption, weighted_labor_supplies, summary_table, ncol=3)
dev.off()

# png(paste0("RA0.1_m0_Lbar_w_change.png"), width=1200, height=900)
# plot_grid(infection, Reff, infect_prob_S, recovereds, total_contacts, hours, labor_S, labor_I, labor_R, ncol=3)
# dev.off()


# png(paste0("RA0.1_m0_planner_zoomedin.png"), width=1200, height=900)
# plot_grid(infection, Reff, infect_prob_S, uninfected, total_contacts, hours, labor_S, labor_I, labor_R, ncol=3)
# dev.off()

# png(paste0("fig2_lines_sketch_all.png"), width=1200, height=900)
# plot_grid(infection, Reff, infect_prob_S, uninfected, total_contacts, hours, labor_S, labor_I, labor_R, ncol=3)
# dev.off()

# png(paste0("fig2_lines_sketch_all_zoomedin.png"), width=1200, height=900)
# plot_grid(infection, Reff, infect_prob_S, uninfected, total_contacts, hours, labor_S, labor_I, labor_R, ncol=3)
# dev.off()

# png(paste0("fig2_lines_sketch_only_econ.png"), width=1200, height=900)
# plot_grid(infection, Reff, infect_prob_S, uninfected, total_contacts, hours, labor_S, labor_I, labor_R, ncol=3)
# dev.off()

# png(paste0("fig2_lines_sketch_all.png"), width=1200, height=900)
# plot_grid(infection, Reff, infect_prob_S, uninfected, total_contacts, hours, labor_S, labor_I, labor_R, ncol=3)
# dev.off()

#####
# Plots for presentation slides
#####


#####
# "supply curve" plots. x-axis is quantity uninfected (susceptible), y-axis is consumption deviations from steady state
#####

# by_type = group_by(all_scenarios, type)

# all_scenarios_summary = summarise(by_type,
# 							neverinfected = min(S),
# 							cost = -min(aggregate_consumption_deviation),
#                                           total_cases = max(total_cases),
#                                           total_deaths = sum(D),
#                                           peak_week_cases = max(weekly_new_cases))

# (MC_within <- ggplot(data=all_scenarios, aes(x=(1-total_cases), group=type, color=type)) + geom_line(aes(y=-aggregate_consumption_deviation), size=1) +
#             theme_bw() + ggtitle("Marginal costs of producing never-infecteds") +
#             xlab("Proportion never infected") + ylab("% C deviation from initial steady state") + 
#             scale_colour_manual(values=Mix))

# (PPF_across <- ggplot(data=all_scenarios_summary, aes(x=neverinfected*100, group=type, color=type)) + 
#             geom_jitter(aes(y=cost), size=4) +
#             theme_bw() + ggtitle("Pandemic possibility space across models:\ntradeoff conclusions you can get from different assumptions
#             	(southeast is optimistic, northwest is pessimistic)") +
#             xlab("Never-infected (% of pop)") + ylab("-% C deviation from initial steady state") + 
#             xlim(-1,100) + ylim(-1,100) +
#             geom_abline(intercept=0,slope=1,linetype="dashed",alpha=0.75,color="darkgray") +
#             scale_colour_manual(values=Mix))

# (model_neverinfecteds = ggplot(data=all_scenarios_summary) + 
#             geom_point(aes(x=neverinfected,y=type, color=type), size=4) +
#             theme_bw() +
#             xlab("Proportion never infected") + ylab("Models") + 
#             scale_colour_manual(values=Mix) + guides(color=FALSE))

# (model_deaths = ggplot(data=all_scenarios_summary) + 
#             geom_point(aes(x=total_deaths,y=type, color=type), size=4) +
#             theme_bw() +
#             xlab("Total deaths") + ylab("Models") + 
#             scale_colour_manual(values=Mix) + guides(color=FALSE))

# (model_peakweek = ggplot(data=all_scenarios_summary) + 
#             geom_point(aes(x=peak_week_cases,y=type, color=type), size=4) +
#             theme_bw() +
#             xlab("Peak week cases") + ylab("Models") + 
#             scale_colour_manual(values=Mix))


# plot_grid(model_neverinfecteds, model_deaths, model_peakweek, ncol=3)

# plot_grid(MC_within,PPF_across)

# png(paste0("../plots/modeling_ppf.png"), width=600, height=400)
# plot_grid(MC_within,PPF_across)
# dev.off()


# png(paste0("fig2_lines_sketch_rebelo_utility_2.png"), width=1200, height=900)
# plot_grid(infection, Reff, cases, uninfected, consumption, hours, labor_S, labor_I, labor_R, model_neverinfecteds, model_deaths, model_peakweek, ncol=3)
# dev.off()


