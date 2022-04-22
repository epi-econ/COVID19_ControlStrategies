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
      geom_contour_filled(aes(x = c_S, y = l_S, z = contacts), show.legend = FALSE) +
      geom_point(aes(x=0,y=0), color="black", fill="white", size=5, pch=21) + 
      geom_point(aes(x=58000/365,y=0), color="black", fill="white", size=5, pch=21) + 
      geom_point(aes(x=0,y=0.3333*24), color="black", fill="white", size=5, pch=21) + 
      # geom_label(aes(x=0,y=0,label="Unavoidable"),hjust=-0.1, vjust=1.25) +
      # geom_label(aes(x=58000/365,y=0,label="Consumption"),hjust=1.08, vjust=1.25) +
      # geom_label(aes(x=0,y=0.3333*24,label="Labor"),hjust=-0.15, vjust=-0.3) +
      labs(y="", x="Consumption by susceptible ($/day)", title="Linear contact function", fill="") +
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
      geom_contour_filled(aes(x = c_S, y = l_S, z = contacts)) +
      geom_point(aes(x=0,y=0), color="black", fill="white", size=5, pch=21) + 
      geom_point(aes(x=58000/365,y=0), color="black", fill="white", size=5, pch=21) + 
      geom_point(aes(x=0,y=0.3333*24), color="black", fill="white", size=5, pch=21) + 
      # geom_label(aes(x=0,y=0,label="Unavoidable"),hjust=-0.1, vjust=1.25) +
      # geom_label(aes(x=58000/365,y=0,label="Consumption"),hjust=1.08, vjust=1.25) +
      # geom_label(aes(x=0,y=0.3333*24,label="Labor"),hjust=-0.15, vjust=-0.3) +
      labs(x="", y="", title="Concave contact function", fill="Contacts/day") +
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
      geom_contour_filled(aes(x = c_S, y = l_S, z = contacts), show.legend = FALSE) +
      geom_point(aes(x=0,y=0), color="black", fill="white", size=5, pch=21) + 
      geom_point(aes(x=58000/365,y=0), color="black", fill="white", size=5, pch=21) + 
      geom_point(aes(x=0,y=0.3333*24), color="black", fill="white", size=5, pch=21) + 
      # geom_label(aes(x=0,y=0,label="Unavoidable"),hjust=-0.1, vjust=1.25) +
      # geom_label(aes(x=58000/365,y=0,label="Consumption"),hjust=1.08, vjust=1.25) +
      # geom_label(aes(x=0,y=0.3333*24,label="Labor"),hjust=-0.15, vjust=-0.3) +
      labs(title="Convex contact function", x="", y="Labor supply by susceptible (hrs/day)", fill="") +
      theme_classic()
gg_convex_contact_surface

contact_surface_plots <- gg_convex_contact_surface | gg_linear_contact_surface | gg_concave_contact_surface
contact_surface_plots

write_csv(convex_contact_surface_data, path="../../../Figures_data_FINAL/figure_si_convex_surface_data.csv")
write_csv(linear_contact_surface_data, path="../../../Figures_data_FINAL/figure_si_linear_surface_data.csv")
write_csv(concave_contact_surface_data, path="../../../Figures_data_FINAL/figure_si_concave_surface_data.csv")

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
# Make plot
#####


png(paste0("../../Results/Figure_SI_contact_function/contact_surface_plot.png"), width=850, height=350)
contact_surface_plots
dev.off()

### Mockup plot with contact surfaces and asymptomatic share

fig4_mockup_2 <- ((asymptomatics_and_productivity / gg_convex_contact_surface / gg_linear_contact_surface / gg_concave_contact_surface) | (recession_ratio_cl | cases_ratio_cl) / (recession_ratio_au | cases_ratio_au) / bar_summary) + plot_annotation(tag_levels = "a")

png(paste0("../../fig4_sketch_with_function_shapes.png"), width=1400, height=1000)
fig4_mockup_2
dev.off()



