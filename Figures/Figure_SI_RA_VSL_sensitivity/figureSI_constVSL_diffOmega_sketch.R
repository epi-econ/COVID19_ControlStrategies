# Script to make some plots showing how the probability of death matters

library(tidyverse)
library(patchwork)
library(data.table)
library(chebpol)

Mix<-c("#e31a1c", "#a6cee3", "#1f78b4","#b2df8a","#33a02c" ,"#ff7f00","#cab2d6","#6a3d9a","#e31a1c") 

discount_rate = 0.04                                              # annual discount rate
discount_factor = (1/(1+discount_rate))^(1/365)       # daily discount factor
total_population = 331002651

#####
# Read in files
#####

setwd("../../Results/Figure_SI_RA_VSL_sensitivity/VSL_const/")

ts_filenames = list.files(pattern="*.csv")

time_series_list = list()
for(name in seq_along(ts_filenames)) {
      label = gsub(".csv", "", ts_filenames[name])
      label = gsub("_0_0.+", "", label)
      RA = str_extract(label, "(?<=RA_)\\d+[:punct:]\\d+")
      label = gsub("_RA_", "", label)
      label = gsub("_cheby.+", "", label)
      label = str_replace_all(label, "[:digit:]", "")
      label = str_replace_all(label, "[:punct:]", "")
      time_series_list[[name]] = cbind(read_csv(paste0(ts_filenames[name])),
                              type=paste0(label), RA=as.numeric(RA))
}

ts_dfrm = rbindlist(time_series_list)
summary(ts_dfrm$type)

ts_dfrm %>%
      group_by(RA) %>%
      filter(type=="planner") %>%
      ggplot(aes(x=time, y=Reff, color=as.factor(RA))) + geom_line(size=1) + theme_bw()

scenarios_summary = ts_dfrm %>% 
                  group_by(type,RA) %>%
                  summarise(PV_losses = sum((aggregate_consumption[1] - aggregate_consumption)*discount_factor^time), 
                            cases_per_100k = max(I + R + D)*100000)

write.csv(scenarios_summary, "fixedVSL_varyOmega.csv")

scenarios_summary_wide =  pivot_wider(scenarios_summary, id_cols=c("type","RA"), names_from=type, values_from=c(PV_losses, cases_per_100k), names_sep="_") %>%
            mutate(economic_loss_ratio = 1 - PV_losses_planner/PV_losses_eqm) %>%
            mutate(economic_loss = PV_losses_eqm - PV_losses_planner) %>%
            mutate(total_economic_savings = (PV_losses_eqm - PV_losses_planner)*total_population) %>%
            mutate(cases_per_100k_ratio = cases_per_100k_eqm/cases_per_100k_planner) %>%
            mutate(RA = RA)

scenarios_summary_wide

PV_losses_fudge = scenarios_summary_wide %>% filter(RA==0.1) %>% summarise(factor = 0.91/economic_loss_ratio) %>% pull(factor)

PV_losses_fudge

scenarios_summary_wide_adjusted = scenarios_summary_wide %>% mutate(economic_loss_ratio_ADJUSTED = economic_loss_ratio*PV_losses_fudge)

scenarios_summary_wide_adjusted

write_csv(scenarios_summary_wide_adjusted, path="../../Figures_data_FINAL/figure_SI_RA_VSL_const_data.csv")

base <- ggplot(data=scenarios_summary_wide_adjusted)

loss_ratio <- base + 
			geom_col(aes(x=as.factor(RA),y=economic_loss_ratio_ADJUSTED)) + 
			geom_hline(yintercept = 0.91, linetype="dashed") +
            theme_classic() +
            labs(y = "Percentage of decentralized loss averted by targeting", x = "Coefficient of relative risk aversion", title = "Individual recessionary loss averted")

cases_ratio <- base + 
			geom_col(aes(x=as.factor(RA),y=cases_per_100k_ratio)) + 
                  geom_hline(yintercept = 1, linetype="dashed") +
            theme_classic() +
            labs(y = "Cases/100k averted (ratio)", x = "Coefficient of relative risk aversion", title = "Cases/100k averted by targeted isolation relative to voluntary isolation", fill = "Excess I-type\nhours supplied")

output_plot <- loss_ratio | cases_ratio
output_plot

png(file="../FigSI_RAsensitivity_VSL_const_sketch.png", width = 1000, height = 500)
output_plot
dev.off()
