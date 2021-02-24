# Script to make some plots for comparing model features. Not meant to be paper-ready, just for internal use

library(ggpubr)
library(tidyverse)
library(patchwork)
library(data.table)
library(janitor)
library(scales)

source("functions.R")

#Mix<-c("#e31a1c", "#a6cee3", "#1f78b4","#b2df8a","#33a02c" ,"#ff7f00","#cab2d6","#6a3d9a","#e31a1c") 
Mix<-c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99") 

parms <- read.csv("../../Results/Figure2_Models_epi/Data_For_plots/value_functions/calibrated_parameters.csv")
parms$kappa_c <- 1
parms$kappa_l <- 1

discount_rate = 0.04                                              # annual discount rate
discount_factor = (1/(1+discount_rate))^(1/365)       # daily discount factor

setwd("../../Results/Figure_SI_blanket_sensitivity/data")

filenames = list.files(pattern="*.csv")

list_of_files = list()
for(name in seq_along(filenames)) {
      label = gsub(".csv", "", filenames[name])
      list_of_files[[name]] = fread(paste0(filenames[name]))
      if(!("ensemble_id" %in% colnames(list_of_files[[name]]))) {
            list_of_files[[name]] = cbind(list_of_files[[name]], ensemble_id = paste0(label))
      }
}

policy_paramset = read.csv("../ensemble_policy_parameters.csv")
policy_paramset$ensemble_id = as.factor(policy_paramset$ensemble_id)

all_scenarios = rbindlist(list_of_files)
all_scenarios = all_scenarios %>% 
                  mutate(weighted_labor_S = S*labor_S) %>%
                  mutate(weighted_labor_I = I*labor_I) %>%
                  mutate(weighted_labor_R = R*labor_R) %>% 
                  mutate(total_contacts = consumption_contacts + labor_contacts + other_contacts) %>%
                  mutate(ensemble_id = recode(ensemble_id, eqm = "decentralized", planner = "coordinated")) %>%
                  mutate(SI_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_I, l_S = labor_S, l_I = labor_I) )) %>%
                  mutate(SR_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_R, l_S = labor_S, l_I = labor_R) )) %>%
                  mutate(RI_contacts = contacts(parms, choices = data.frame(c_S = consumption_R, c_I = consumption_I, l_S = labor_R, l_I = labor_I) )) %>%
                  mutate(SS_contacts = contacts(parms, choices = data.frame(c_S = consumption_S, c_I = consumption_S, l_S = labor_S, l_I = labor_S) )) %>%
                  mutate(RR_contacts = contacts(parms, choices = data.frame(c_S = consumption_R, c_I = consumption_R, l_S = labor_R, l_I = labor_R) )) %>%
                  mutate(prob_contact_I = (SI_contacts + RI_contacts)/(SI_contacts + RI_contacts + SR_contacts + SS_contacts + RR_contacts)) %>%
                  mutate(prob_contact_I_weighted = (SI_contacts*S*I + RI_contacts*R*I)/(SI_contacts*S*I + RI_contacts*R*I + SR_contacts*S*R + SS_contacts*S*S + RR_contacts*R*R))

scenarios_summary = all_scenarios %>% 
                  group_by(ensemble_id) %>%
                  summarise(PV_losses = sum(((58000/365) - aggregate_consumption)*discount_factor^time), 
                            cases_per_100k = max(I + R + D)*100000) %>%
                  left_join(policy_paramset) %>%
                  mutate(duration = end_date - start_date) %>%
                  mutate(intensity = duration*(1-labor_supply_level))
scenarios_summary

q_cases <- quantile(scenarios_summary$cases_per_100k)
q_loss <- quantile(scenarios_summary$PV_losses)

scenarios_summary$q_cases <- cut(scenarios_summary$cases_per_100k, q_cases)
levels(scenarios_summary$q_cases) <- c("Q1", "Q2", "Q3", "Q4")
scenarios_summary

scenarios_summary$q_losses <- cut(scenarios_summary$PV_losses, q_losses)
levels(scenarios_summary$q_losses) <- c("Q1", "Q2", "Q3", "Q4")
scenarios_summary

cases_25th_percentile <- which.max(scenarios_summary$cases_per_100k[scenarios_summary$q_cases=="Q1"])
cases_75th_percentile <- which.max(scenarios_summary$cases_per_100k[scenarios_summary$q_cases=="Q3"])

# scenarios_summary_annotated = left_join(scenarios_summary, policy_paramset, by=c("ensemble_id"))

scenarios_summary[which(scenarios_summary$cases_per_100k==as.numeric(quantile(scenarios_summary$cases_per_100k, probs=c(0.25)))),]


write.csv(scenarios_summary, file="../ensemble_summary.csv")

## Make plots

line_data = all_scenarios %>% filter(ensemble_id != "decentralized" , ensemble_id !="coordinated")

lines_base = ggplot(data=line_data, aes(x=time))

infection = ggplot() + 
            geom_line(data = line_data, aes(x = time, y=I*100, group=ensemble_id), alpha = 0.05, size = 0.75) +
            geom_line(data =(all_scenarios %>% filter(ensemble_id == "decentralized")), aes(x = time, y = I*100), size=1, color="firebrick3") +
            geom_line(data =(all_scenarios %>% filter(ensemble_id == "coordinated")), aes(x = time, y = I*100), linetype="dashed", size=1.5, color="dodgerblue3") +
            theme_bw() + ggtitle("Infecteds") +
            xlab("Day") + ylab("Proportion of population (%)") + 
            guides(linetype=FALSE, size=FALSE)

deaths = ggplot() + 
            geom_line(data = line_data, aes(x = time, y=D*100, group=ensemble_id), alpha = 0.05, size = 0.75) +
            geom_line(data =(all_scenarios %>% filter(ensemble_id == "decentralized")), aes(x = time, y = D*100), size=1, color="firebrick3") +
            geom_line(data =(all_scenarios %>% filter(ensemble_id == "coordinated")), aes(x = time, y = D*100), linetype="dashed", size=1.5, color="dodgerblue3") +
            theme_bw() + ggtitle("Cumulative deaths") +
            xlab("Day") + ylab("Proportion of population (%)") + 
            guides(linetype=FALSE, size=FALSE)

recovereds = ggplot() + 
            geom_line(data = line_data, aes(x = time, y=R*100, group=ensemble_id), alpha = 0.05, size = 0.75) +
            geom_line(data =(all_scenarios %>% filter(ensemble_id == "decentralized")), aes(x = time, y = R*100), size=1, color="firebrick3") +
            geom_line(data =(all_scenarios %>% filter(ensemble_id == "coordinated")), aes(x = time, y = R*100), linetype="dashed", size=1.5, color="dodgerblue3") +
            theme_bw() + ggtitle("Recovereds") +
            xlab("Day") + ylab("Proportion of population (%)") + 
            guides(linetype=FALSE, size=FALSE)

R_effective = ggplot() + 
            geom_line(data = line_data, aes(x = time, y=Reff*100, group=ensemble_id), alpha = 0.05, size = 0.75) +
            geom_line(data =(all_scenarios %>% filter(ensemble_id == "decentralized")), aes(x = time, y = Reff*100), size=1, color="firebrick3") +
            geom_line(data =(all_scenarios %>% filter(ensemble_id == "coordinated")), aes(x = time, y = Reff*100), linetype="dashed", size=1.5, color="dodgerblue3") +
            theme_bw() + ggtitle("R_effective") +
            xlab("Day") + ylab("R_t") + 
            guides(linetype=FALSE, size=FALSE)

susceptibles = ggplot() + 
            geom_line(data = line_data, aes(x = time, y=S*100, group=ensemble_id), alpha = 0.05, size = 0.75) +
            geom_line(data =(all_scenarios %>% filter(ensemble_id == "decentralized")), aes(x = time, y = S*100), size=1, color="firebrick3") +
            geom_line(data =(all_scenarios %>% filter(ensemble_id == "coordinated")), aes(x = time, y = S*100), linetype="dashed", size=1.5, color="dodgerblue3") +
            theme_bw() + ggtitle("Susceptibles") +
            xlab("Day") + ylab("Proportion of population (%)") + 
            guides(linetype=FALSE, size=FALSE)

consumption = ggplot() + 
            geom_line(data = line_data, aes(x = time, y=aggregate_consumption_deviation, group=ensemble_id), alpha = 0.05, size = 0.75) +
            geom_line(data =(all_scenarios %>% filter(ensemble_id == "decentralized")), aes(x = time, y = aggregate_consumption_deviation), size=1, color="firebrick3") +
            geom_line(data =(all_scenarios %>% filter(ensemble_id == "coordinated")), aes(x = time, y = aggregate_consumption_deviation), linetype="dashed", size=1.5, color="dodgerblue3") +
            theme_bw() + ggtitle("GDP deviation") +
            xlab("Day") + ylab("Percentage of initial level") + 
            guides(linetype=FALSE, size=FALSE)

hist_data = scenarios_summary %>% filter(ensemble_id != "decentralized" , ensemble_id !="coordinated")
n_ensembles = nrow(hist_data)

epi_hist = ggplot(data = hist_data, aes(x=cases_per_100k)) + 
            geom_vline(xintercept= as.numeric(scenarios_summary %>% filter(ensemble_id=="decentralized") %>% summarise(cases_per_100k) ), size=1.5, color="firebrick3") +
            geom_vline(xintercept=as.numeric(scenarios_summary %>% filter(ensemble_id=="coordinated") %>% summarise(cases_per_100k) ), size=1.5, linetype="dashed", color="dodgerblue3" ) + 
            geom_histogram(bins = 100) +
            #xlim(min(scenarios_summary$cases_per_100k),max(scenarios_summary$cases_per_100k)) +
            theme_bw() + ggtitle("Total cases per 100,000") +
            xlab("Cases/100k")

econ_hist = ggplot(data = hist_data, aes(x=PV_losses)) + 
            geom_vline(xintercept= as.numeric(scenarios_summary %>% filter(ensemble_id=="decentralized") %>% summarise(PV_losses) ), size=1.5, color="firebrick3") +
            geom_vline(xintercept=as.numeric(scenarios_summary %>% filter(ensemble_id=="coordinated") %>% summarise(PV_losses) ), size=1.5, linetype="dashed" , color="dodgerblue3") +
            geom_histogram(bins = 100) + 
            xlim(min(scenarios_summary$PV_losses),max(scenarios_summary$PV_losses)) +
            theme_bw() + ggtitle("PV of average individual losses") +
            xlab("PV $ lost for average individual")

policy_frontier = ggplot() + 
            geom_point(data = hist_data, aes(x=cases_per_100k, y=PV_losses), alpha=0.075) +
            geom_point(data = (scenarios_summary %>% filter(ensemble_id=="decentralized")), aes(x=cases_per_100k, y=PV_losses), color="firebrick3", size = 2) +
            geom_point(data = (scenarios_summary %>% filter(ensemble_id=="coordinated")), aes(x=cases_per_100k, y=PV_losses), color="dodgerblue3", size = 2) +
            theme_bw() + ggtitle("Policy efficiency space (further southwest is better)") +
            xlab("Cases/100k") + ylab("PV $ lost for average individual") + coord_flip()

lockdown_frontier = ggplot() + 
            geom_point(data = hist_data, aes(x=cases_per_100k, y=PV_losses, color= intensity)) +
            theme_bw() + ggtitle("Policy efficiency space (blanket lockdowns only)") +
            xlab("Cases/100k") + ylab("PV $ lost for average individual") + coord_flip()

ensemble_plot = (infection | consumption) / ( R_effective | susceptibles | recovereds | deaths) / (epi_hist | econ_hist) / (policy_frontier | lockdown_frontier) + 
      plot_layout(heights = c(2,2,1,2)) + 
      plot_annotation(
        title = 'The relative inefficacy of blanket lockdowns',
        subtitle = paste0('Ensemble of ', n_ensembles, ' random blanket policies. Red is decentralized, blue is coordinated, black is general lockdown.'),
        caption = 'Optimally-targeted isolation policy produces the smallest recession and among the best disease outcomes.'
      )

ensemble_plot

png(paste0("../ensemble_plot_full.png"), width=1500, height=900)
ensemble_plot
dev.off()

ensemble_plot_main <- (policy_frontier | lockdown_frontier) + plot_annotation(tag_levels = 'a')

png(paste0("../lockdown_main.png"), width=1000, height=500)
ensemble_plot_main
dev.off()

ensemble_plot_second <- (consumption | infection) / (econ_hist | epi_hist) + plot_annotation(tag_levels = 'a')

png(paste0("../lockdown_second.png"), width=1000, height=500)
ensemble_plot_second
dev.off()


ensemble_plot_small = (consumption | infection) / (econ_hist | epi_hist) / (policy_frontier | lockdown_frontier) + plot_annotation(tag_levels = 'a')

ensemble_plot_small

png(paste0("../ensemble_plot_small.png"), width=1500, height=900)
ensemble_plot_small
dev.off()


