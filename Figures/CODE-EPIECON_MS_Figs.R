#######################################################
##############  Manuescript figures Figures ####################
#######################################################




#######################################################
##############  Conceptual figure ####################
#######################################################

# FIGURE 1
#Conceptual figure made in illustrator




#######################################################
##############  Dynamics ####################
#######################################################



# FIGURE 2



#a) Epi evolution
#libraries
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(gridExtra)
library(hrbrthemes)

#read in data
Fig2<-read.csv("/~/Fig2/figure_2_lineplot_data.csv", header=TRUE,sep=",")

#plot main


#cols<-c("#7fc97f","#bf5b17", "#386cb0")

cols1 = c("#d95f02","#1b9e77","#7570b3")
 p1<-ggplot(Fig2, aes(x=time, y=I, fill=type, color=type)) + geom_area(position = "identity", alpha=0.5, size=0.6) + scale_fill_manual(values=cols, name = "Policy") + scale_color_manual(values=cols, name = "Policy")+
theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),legend.position="none")+ylab("Proportion of cases")+xlab("Day")+scale_x_continuous(breaks=seq(0, 600, 100))




#b) losses barplot

loss <- data.frame(pol=c("A", "B", "C"),
                value=c(21119.23, 12078.67,1082.67))
                
p2<-ggplot(data=loss, aes(x=pol, y=value, fill=pol))+
  geom_bar(stat="identity")+scale_fill_manual(values=c("#d95f02","#1b9e77","#7570b3")) +scale_x_discrete(breaks=c("A","B","C"),
        labels=c("Blanket \n Lockdown","Voluntary \n Isolation","Targeted \n Isolation"))+theme_minimal()+ylab("Avg. Individual losses \n ($/person)")+xlab("") + theme_ipsum(axis_title_size = 24, base_size = 20)+theme(legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+coord_cartesian(ylim=c(0,24000))+scale_y_continuous(breaks=seq(0,24000, 8000))              

#c) summary table
## This part is not coded in R -- assembled in illustrator 


## See Fig parts final figure assembled in illustrator



#######################################################################################################################        
#######################################################################################################################
###################################   MECHANISMS  ###########################################################  
####################################################################################################################### 


# FIGURE  3

library("wesanderson")

cols1 = c("#1b9e77","#7570b3") #"#d95f02"


 Zis = c("#F21A00","#EBCC2A",,"#3B9AB2") #"#78B7C5",,"#E1AF00"
Fig3a<-read.csv("/~/Fig3/Fig3a.csv", header=TRUE,sep=",")
#3a
###If done with three colors by sidease status
p<-ggplot(Fig3a, aes(x=time, y=value, colour=Labor, linetype=type))+
    geom_line(size=0.8,alpha=0.7, position = position_dodge(width = 0.8)) +
scale_colour_manual(values=Zis, name = "Disease Status") +
    theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major =element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),legend.position="none")+ylab("Normalized person-hours")+xlab("Day")+scale_x_continuous(breaks=seq(0, 600, 100)) +scale_linetype_manual(values = c( 1,3), name = "Policy") +guides(linetype=FALSE)      



# #3a

# p<-ggplot(Fig3a, aes(x=time, y=value, colour=type, linetype=Labor))+
    # geom_line(size=1, alpha=0.8) +
     # scale_colour_manual(values=cols1, name = "Policy") +
    # theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major =element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),legend.position="none")+ylab("Normalized person-hours (sqrt)")+xlab("Day")+scale_x_continuous(breaks=seq(0, 600, 100))+coord_cartesian(ylim=c(0,8))+scale_y_continuous(breaks=seq(0, 8, 1)) +scale_linetype_manual(values = c( 1,3, 2),name = "Disease status")+guides(linetype=FALSE) 



#3b
Fig3<-read.csv("/~Fig3/figure_2_lineplot_data.csv", header=TRUE,sep=",")
p1<-ggplot(Fig3, aes(x=time, y=prob_contact_I_weighted, colour=type)) + 
    geom_line(size=0.8, alpha=0.9) +  scale_colour_manual(values=cols1, name = "Policy") +
    theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major =element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),legend.position="none")+ylab("Probability S/R contacts I")+xlab("Day")+coord_cartesian(ylim=c(0,0.04))+scale_x_continuous(breaks=seq(0, 600, 100))+scale_y_continuous(breaks=seq(0, 0.04, 0.01))


#3c SI contacts

p2<-ggplot(Fig3, aes(x=time, y=total_contacts, colour=type,linetype = type)) + 
 geom_line(size=1) +
    scale_colour_manual(values=cols1, name = "Policy") +
    theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),legend.position="none")+ylab("Mean daily contacts S-I")+xlab("Day")+scale_x_continuous(breaks=seq(0, 600, 100))
        
#3d Avg contacts by activity


Fig3de<-read.csv("/~/Fig3/Fig3de.csv", header=TRUE,sep=",")

p3<-ggplot(Fig3de, aes(x=time, y=contacts, colour=type, linetype=policy))+
    geom_line(size=1, alpha=0.8) +
     scale_colour_manual(values=cols, name = "Activity") +
    theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major =element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),,legend.position="none")+ylab("Average contacts (activity type)")+xlab("Day")+scale_x_continuous(breaks=seq(0, 600, 100))+coord_cartesian(ylim=c(0,8))+scale_y_continuous(breaks=seq(0, 8, 1)) +scale_linetype_manual(values = c( 1,3, 2),name = "Policy")+guides(linetype=FALSE) 

#3e Cases by activity


p4<-ggplot(Fig3de, aes(x=time, y=cases, colour=type, linetype=policy))+
    geom_line(size=1, alpha=0.8) +
     scale_colour_manual(values=cols, name = "Activity") +
    theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major =element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),legend.position="none")+ylab("Prevalence (activity type)")+xlab("Day")+scale_x_continuous(breaks=seq(0, 600, 100))+coord_cartesian(ylim=c(0,0.0025))+scale_y_continuous(breaks=seq(0, 0.0025, 0.001)) +scale_linetype_manual(values = c( 1,3, 2),name = "Policy")+guides(linetype=FALSE) 
    
#######################################################################################################################        
#######################################################################################################################
###################################   FRICTIONS  ###########################################################  
####################################################################################################################### 

# FIGURE 4

#info lag figure
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(gridExtra)
library(hrbrthemes)

Figinfoleft<-read.csv("/~/FigureX_infolag/linear_behavior/Worst-case_dynamics.csv", header=TRUE,sep=",") 
Figinforight<-read.csv("/~/FigureX_infolag/linear_behavior/Worst-case_summary_stats.csv", header=TRUE,sep=",") 

#Row 1 Lowest test quality
cols1 = c("#d95f02","#7570b3","#1b9e77")

## using area
#p1<-ggplot(Figinfoleft, aes(x=time, y=I*100, fill=type, color=type)) + geom_area(position = "identity", alpha=0.5, size=0.6) + scale_fill_manual(values=cols1, name = "Policy") + scale_color_manual(values=cols1, name = "Policy")+theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),legend.position="none")+ylab("Proportion of cases")+xlab("Day")+scale_x_continuous(breaks=seq(0, 250, 50))

# LIne plots
p1<-ggplot(Figinfoleft, aes(x=time, y=I*100, color=type)) +  geom_line(size=1, alpha=0.9) +  scale_colour_manual(values=cols1, name = "Policy")+theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),legend.position="none")+ylab("% current infections")+xlab("Day")+scale_x_continuous(breaks=seq(0, 250, 50))+coord_cartesian(ylim=c(0,10))


p2<-ggplot(Figinfoleft, aes(x=time, y=aggregate_consumption_deviation, fill=type, color=type)) + geom_line(size=1, alpha=0.9)+ scale_color_manual(values=cols1, name = "Policy")+theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),legend.position="none")+ylab("% Deviation")+xlab("Day")+scale_x_continuous(breaks=seq(0, 250,50))+coord_cartesian(ylim=c(-80, 0))



#b) losses barplot

loss <- data.frame(pol=c("A", "B", "C"),
                value=c(8846, 5485,4814))
                
p3<-ggplot(data=loss, aes(x=pol, y=value, fill=pol))+
  geom_bar(stat="identity")+scale_fill_manual(values=c("#d95f02","#1b9e77","#7570b3")) +scale_x_discrete(breaks=c("A","B","C"),
        labels=c("Blanket \n Lockdown","Voluntary \n Isolation","Targeted \n Isolation"))+theme_minimal()+ylab("Avg.loss ($/person)")+xlab("") + theme_ipsum(axis_title_size = 24, base_size = 20)+theme(legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+coord_cartesian(ylim=c(0,14000))+scale_y_continuous(breaks=seq(0,14000, 2000))           


cases <- data.frame(pol=c("A", "B", "C"),
                value=c(63110, 63170,63109))
p4<-ggplot(data=cases, aes(x=pol, y=value, color=pol))+
  geom_point(size=4)+scale_color_manual(values=c("#d95f02","#1b9e77","#7570b3")) +scale_x_discrete(breaks=c("A","B","C"),
        labels=c("Blanket \n Lockdown","Voluntary \n Isolation","Targeted \n Isolation"))+theme_minimal()+ylab("Prevalence / 100,000")+xlab("") + theme_ipsum(axis_title_size = 24, base_size = 20)+theme(legend.position = "none",)+coord_cartesian(ylim=c(60000,64000))+scale_y_continuous(breaks=seq(60000,64000, 1000)) 

worst-case_dynamics.csv  2 #on the right 
worst-case_summary_stats.csv #2 on the left


#row 2 improved test quality


Figinfoleft2<-read.csv("/~/FigureX_infolag/linear_behavior/Improvingtestquality_dynamics.csv", header=TRUE,sep=",") 
Figinforight2<-read.csv("/~/FigureX_infolag/linear_behavior/Worst-case_summary_stats.csv", header=TRUE,sep=",") 

# LIne plots
p1<-ggplot(Figinfoleft2, aes(x=time, y=I*100, color=type)) +  geom_line(size=1, alpha=0.9) +  scale_colour_manual(values=cols1, name = "Policy")+theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),legend.position="none")+ylab("% current infections")+xlab("Day")+scale_x_continuous(breaks=seq(0, 250, 50))+coord_cartesian(ylim=c(0,10))


p2<-ggplot(Figinfoleft2, aes(x=time, y=aggregate_consumption_deviation, fill=type, color=type)) + geom_line(size=1, alpha=0.9)+ scale_color_manual(values=cols1, name = "Policy")+theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),legend.position="none")+ylab("% Deviation")+xlab("Day")+scale_x_continuous(breaks=seq(0, 250,50))+coord_cartesian(ylim=c(-80, 0))



#b) losses barplot

loss2 <- data.frame(pol=c("A", "B", "C"),
                 value=c(12502, 9891,1408))
                
p3<-ggplot(data=loss2, aes(x=pol, y=value, fill=pol))+
  geom_bar(stat="identity")+scale_fill_manual(values=c("#d95f02","#1b9e77","#7570b3")) +scale_x_discrete(breaks=c("A","B","C"),
        labels=c("Blanket \n Lockdown","Voluntary \n Isolation","Targeted \n Isolation"))+theme_minimal()+ylab("Avg.loss ($/person)")+xlab("") + theme_ipsum(axis_title_size = 24, base_size = 20)+theme(legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+coord_cartesian(ylim=c(0,14000))+scale_y_continuous(breaks=seq(0,14000, 2000))              


cases <- data.frame(pol=c("A", "B", "C"),
                value=c(61480, 61545,61347))
p4<-ggplot(data=cases, aes(x=pol, y=value, color=pol))+
  geom_point(size=4)+scale_color_manual(values=c("#d95f02","#1b9e77","#7570b3")) +scale_x_discrete(breaks=c("A","B","C"),
        labels=c("Blanket \n Lockdown","Voluntary \n Isolation","Targeted \n Isolation"))+theme_minimal()+ylab("Prevalence / 100,000")+xlab("") + theme_ipsum(axis_title_size = 24, base_size = 20)+theme(legend.position = "none",)+coord_cartesian(ylim=c(60000,64000))+scale_y_continuous(breaks=seq(60000,64000, 1000)) 

worst-case_dynamics.csv  2 #on the right 
worst-case_summary_stats.csv #2 on the left



# row 3

Figinfoleft3<-read.csv("/~/FigureX_infolag/linear_behavior/Improvingtestqualityandtimeliness_dynamics.csv", header=TRUE,sep=",") 


# LIne plots
p1<-ggplot(Figinfoleft3, aes(x=time, y=I*100, color=type)) +  geom_line(size=1, alpha=0.9) +  scale_colour_manual(values=cols1, name = "Policy")+theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),legend.position="none")+ylab("% current infections")+xlab("Day")+scale_x_continuous(breaks=seq(0, 250, 50))+coord_cartesian(ylim=c(0,10))


p2<-ggplot(Figinfoleft3, aes(x=time, y=aggregate_consumption_deviation, fill=type, color=type)) + geom_line(size=1, alpha=0.9)+ scale_color_manual(values=cols1, name = "Policy")+theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),legend.position="none")+ylab("% Deviation")+xlab("Day")+scale_x_continuous(breaks=seq(0, 250,50))+coord_cartesian(ylim=c(-80, 0))






#b) losses barplot

loss3 <- data.frame(pol=c("A", "B", "C"),
                 value=c(12434, 9615,1409))
                
p3<-ggplot(data=loss3, aes(x=pol, y=value, fill=pol))+
  geom_bar(stat="identity")+scale_fill_manual(values=c("#d95f02","#1b9e77","#7570b3")) +scale_x_discrete(breaks=c("A","B","C"),
        labels=c("Blanket \n Lockdown","Voluntary \n Isolation","Targeted \n Isolation"))+theme_minimal()+ylab("Avg.loss ($/person)")+xlab("") + theme_ipsum(axis_title_size = 24, base_size = 20)+theme(legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+coord_cartesian(ylim=c(0,14000))+scale_y_continuous(breaks=seq(0,14000, 2000))              


cases3 <- data.frame(pol=c("A", "B", "C"),
                value=c(61603, 61663,61466))
p4<-ggplot(data=cases3, aes(x=pol, y=value, color=pol))+
  geom_point(size=6)+scale_color_manual(values=c("#d95f02","#1b9e77","#7570b3")) +scale_x_discrete(breaks=c("A","B","C"),
        labels=c("Blanket \n Lockdown","Voluntary \n Isolation","Targeted \n Isolation"))+theme_minimal()+ylab("Prevalence / 100,000")+xlab("") + theme_ipsum(axis_title_size = 24, base_size = 20)+theme(legend.position = "none",)+coord_cartesian(ylim=c(60000,64000))+scale_y_continuous(breaks=seq(60000,64000, 1000)) 



#row 4 
library(RColorBrewer)

Figinfoleft4<-read.csv("/~FigureX_infolag/linear_behavior/infolag_scenarios_summary_ratios.csv", header=TRUE,sep=",") 
#infolag_secanrios_summary_ratios.csv    #bottom
legend_title <- "Scenarios"
p5<-ggplot(data=Figinfoleft4, aes(x=perc_change_loss, y=perc_change_case, fill=type))+geom_abline(slope=1, linetype="dotted", color="grey60")+geom_point(shape = 21, colour = "grey20", size = 5, stroke = 0.7)+scale_fill_brewer(palette = "GnBu",labels=c("Testing lag + poor quality test","Testing lag + improving test quality","Decrease testing lag + improving test quality","Frictionless baseline"), legend_title)+theme_minimal()+theme_ipsum(axis_title_size = 14, base_size = 14)+theme(legend.text = element_text(size=14))+coord_cartesian(ylim=c(0,1),xlim=c(0,1))+scale_y_continuous(breaks=seq(0,1, 0.25)) +ylab("Max disease control gains achieved (%)")+xlab("Max economic gains achieved (%)")

## assemblaged in illustrator

#################################################################
# FIGURE 5

# compliance lag figure
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(gridExtra)
library(hrbrthemes)

Figinfoleft<-read.csv("/~/FigureX_compliance/linear_behavior/Zerocompliance-perfectinformation_dynamics.csv", header=TRUE,sep=",") 


#Row 1 Lowest test quality
cols1 = c("#d95f02","#1b9e77","#7570b3")



# LIne plots
p1<-ggplot(Figinfoleft, aes(x=time, y=I*100, color=type)) +  geom_line(size=1, alpha=0.9) +  scale_colour_manual(values=cols1, name = "Policy")+theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),legend.position="none")+ylab("% current infections")+xlab("Day")+scale_x_continuous(breaks=seq(0, 250, 50))+coord_cartesian(ylim=c(0,10))


p2<-ggplot(Figinfoleft, aes(x=time, y=aggregate_consumption_deviation, fill=type, color=type)) + geom_line(size=1, alpha=0.9)+ scale_color_manual(values=cols1, name = "Policy")+theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),legend.position="none")+ylab("% Deviation")+xlab("Day")+scale_x_continuous(breaks=seq(0, 250,50))+coord_cartesian(ylim=c(-80, 0))





#b) losses barplot

loss <- data.frame(pol=c("A", "B", "C"),
                value=c(12033, 12033,12033))
                
p3<-ggplot(data=loss, aes(x=pol, y=value, fill=pol))+
  geom_bar(stat="identity")+scale_fill_manual(values=c("#d95f02","#1b9e77","#7570b3")) +scale_x_discrete(breaks=c("A","B","C"),
        labels=c("Blanket \n Lockdown","Voluntary \n Isolation","Targeted \n Isolation"))+theme_minimal()+ylab("Avg.loss ($/person)")+xlab("") + theme_ipsum(axis_title_size = 24, base_size = 20)+theme(legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+coord_cartesian(ylim=c(0,14000))+scale_y_continuous(breaks=seq(0,14000, 2000))           


cases <- data.frame(pol=c("A", "B", "C"),
                value=c(61792, 61792,61792))
p4<-ggplot(data=cases, aes(x=pol, y=value, color=pol))+
  geom_point(size=4)+scale_color_manual(values=c("#d95f02","#1b9e77","#7570b3")) +scale_x_discrete(breaks=c("A","B","C"),
        labels=c("Blanket \n Lockdown","Voluntary \n Isolation","Targeted \n Isolation"))+theme_minimal()+ylab("Prevalence / 100,000")+xlab("") + theme_ipsum(axis_title_size = 24, base_size = 20)+theme(legend.position = "none",)+coord_cartesian(ylim=c(60000,64000))+scale_y_continuous(breaks=seq(60000,64000, 1000)) 

worst-case_dynamics.csv  2 #on the right 
worst-case_summary_stats.csv #2 on the left


#row 2 improved test quality


Figinfoleft2<-read.csv("/~/FigureX_compliance/linear_behavior/Partialcompliance-improvinginformation_dynamics.csv", header=TRUE,sep=",") 


# LIne plots
p1<-ggplot(Figinfoleft2, aes(x=time, y=I*100, color=type)) +  geom_line(size=1, alpha=0.9) +  scale_colour_manual(values=cols1, name = "Policy")+theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),legend.position="none")+ylab("% current infections")+xlab("Day")+scale_x_continuous(breaks=seq(0, 250, 50))+coord_cartesian(ylim=c(0,10))


p2<-ggplot(Figinfoleft2, aes(x=time, y=aggregate_consumption_deviation, fill=type, color=type)) + geom_line(size=1, alpha=0.9)+ scale_color_manual(values=cols1, name = "Policy")+theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),legend.position="none")+ylab("% Deviation")+xlab("Day")+scale_x_continuous(breaks=seq(0, 250,50))+coord_cartesian(ylim=c(-80, 0))




#b) losses barplot

loss2 <- data.frame(pol=c("A", "B", "C"),
                 value=c(12434, 9615,1409))
                
p3<-ggplot(data=loss2, aes(x=pol, y=value, fill=pol))+
  geom_bar(stat="identity")+scale_fill_manual(values=c("#d95f02","#1b9e77","#7570b3")) +scale_x_discrete(breaks=c("A","B","C"),
        labels=c("Blanket \n Lockdown","Voluntary \n Isolation","Targeted \n Isolation"))+theme_minimal()+ylab("Avg.loss ($/person)")+xlab("") + theme_ipsum(axis_title_size = 24, base_size = 20)+theme(legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+coord_cartesian(ylim=c(0,14000))+scale_y_continuous(breaks=seq(0,14000, 2000))              


cases <- data.frame(pol=c("A", "B", "C"),
                value=c(61603, 61663,61466))
p4<-ggplot(data=cases, aes(x=pol, y=value, color=pol))+
  geom_point(size=4)+scale_color_manual(values=c("#d95f02","#1b9e77","#7570b3")) +scale_x_discrete(breaks=c("A","B","C"),
        labels=c("Blanket \n Lockdown","Voluntary \n Isolation","Targeted \n Isolation"))+theme_minimal()+ylab("Prevalence / 100,000")+xlab("") + theme_ipsum(axis_title_size = 24, base_size = 20)+theme(legend.position = "none",)+coord_cartesian(ylim=c(60000,64000))+scale_y_continuous(breaks=seq(60000,64000, 1000)) 

#worst-case_dynamics.csv  2 #on the right 
#worst-case_summary_stats.csv #2 on the left




 
# row 2

Figinfoleft3<-read.csv("/~/FigureX_compliance/linear_behavior/Partialcompliance-perfectinformation_dynamics.csv", header=TRUE,sep=",") 


# LIne plots
p1<-ggplot(Figinfoleft3, aes(x=time, y=I*100, color=type)) +  geom_line(size=1, alpha=0.9) +  scale_colour_manual(values=cols1, name = "Policy")+theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),legend.position="none")+ylab("% current infections")+xlab("Day")+scale_x_continuous(breaks=seq(0, 250, 50))+coord_cartesian(ylim=c(0,10))


p2<-ggplot(Figinfoleft3, aes(x=time, y=aggregate_consumption_deviation, fill=type, color=type)) + geom_line(size=1, alpha=0.9)+ scale_color_manual(values=cols1, name = "Policy")+theme_ipsum(axis_title_size = 24, base_size = 24)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(size = (1), colour="grey92"),panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),legend.position="none")+ylab("% Deviation")+xlab("Day")+scale_x_continuous(breaks=seq(0, 250,50))+coord_cartesian(ylim=c(-80, 0))





#b) losses barplot

loss3 <- data.frame(pol=c("A", "B", "C"),
                 value=c(12360, 12033,3819))
                
p3<-ggplot(data=loss3, aes(x=pol, y=value, fill=pol))+
  geom_bar(stat="identity")+scale_fill_manual(values=c("#d95f02","#1b9e77","#7570b3")) +scale_x_discrete(breaks=c("A","B","C"),
        labels=c("Blanket \n Lockdown","Voluntary \n Isolation","Targeted \n Isolation"))+theme_minimal()+ylab("Avg.loss ($/person)")+xlab("") + theme_ipsum(axis_title_size = 24, base_size = 20)+theme(legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+coord_cartesian(ylim=c(0,14000))+scale_y_continuous(breaks=seq(0,14000, 2000))              


cases3 <- data.frame(pol=c("A", "B", "C"),
                value=c(61778,61792,61752))
p4<-ggplot(data=cases3, aes(x=pol, y=value, color=pol))+
  geom_point(size=6)+scale_color_manual(values=c("#d95f02","#1b9e77","#7570b3")) +scale_x_discrete(breaks=c("A","B","C"),
        labels=c("Blanket \n Lockdown","Voluntary \n Isolation","Targeted \n Isolation"))+theme_minimal()+ylab("Prevalence / 100,000")+xlab("") + theme_ipsum(axis_title_size = 24, base_size = 20)+theme(legend.position = "none",)+coord_cartesian(ylim=c(60000,64000))+scale_y_continuous(breaks=seq(60000,64000, 1000)) 

worst-case_dynamics.csv  2 #on the right 
worst-case_summary_stats.csv #2 on the left



#row 4
library(RColorBrewer)

Figinfoleft4<-read.csv("/~/FigureX_compliance/linear_behavior/compliance_scenarios_summary_ratios.csv", header=TRUE,sep=",") 
#infolag_secanrios_summary_ratios.csv   

 #bottom

p5<-ggplot(data=Figinfoleft4, aes(x=perc_change_loss, y=perc_change_case, fill=type))+geom_abline(slope=1, linetype="dotted", color="grey60")+geom_point(shape = 21, colour = "grey20", size = 5, stroke = 0.7)+scale_fill_brewer(palette = "GnBu",labels=c("Zero compliance + perfect information","Partial compliance + perfect information","Partial compliance + improving information","Frictionless baseline"),legend_title)+theme_minimal()+theme_ipsum(axis_title_size = 14, base_size = 14)+theme(legend.text = element_text(size=14))+coord_cartesian(ylim=c(0,1),xlim=c(0,1))+scale_y_continuous(breaks=seq(0,1, 0.25)) +ylab("Max disease control gains achieved (%)")+xlab("Max economic gains achieved (%)")
#compliance


## assemblaged in illustrator

#######################################################################################################################        
#######################################################################################################################
###################################   MAPPINGS ###########################################################  
#######################################################################################################################    
    
 #FIGURE 6
 ### mapping
 Fig4map<-read.csv("/~/Fig4/mapping.csv", header=TRUE,sep=",")
 
 ggplot(Fig4map, aes(x=productivity_loss, y= asymptomatic))+geom_line(colour="orangered", size=1)+theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major =element_line(size=0.5, colour="grey92"),panel.grid.minor.y = element_blank())+ylab("Proportion of pre-symptomatic, asymptomatic or mild infections")+xlab("Productivity loss")+scale_x_continuous(breaks=seq(0, 0.8, 0.1), expand = c(0.005,-0.002)) +scale_y_continuous(breaks=seq(0, 1, 0.25))
 
 ##Heatmaps 

## check the data
 Fig4c<-read.csv("/~/Fig4/NEW/4ac.csv", header=TRUE,sep=",") 

 p<-ggplot(Fig4c) +
  geom_contour_filled(aes(y=value,x=asymptomatic_prop, z = economic_loss_ratio_ADJUSTED),size=0.009)+facet_grid(~ratio)+scale_fill_manual(values = wes_palette(14, name = "Zissou1", type = "continuous"), name = "Loss averted (%)")+theme_minimal()+ylab("Ratio of consumption-to-labor contacts") +xlab("Productivity loss (%)") +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_text(color="grey30"),legend.position = "none")+scale_y_continuous(breaks = seq(0,1.25,by = 0.25))+scale_x_continuous(breaks = seq(0,90,by = 10))
    #+geom_point(aes(x=14.45368, y=4.2877), fill="white", colour="#4f418b", shape=21,size = 1.5, stroke =1) 
    
 #4b
 Fig4a<-read.csv("/~/Fig4/figure_4_cl_heatmap_data.csv", header=TRUE,sep=",")   
 
p<-ggplot(Fig4a) +
  geom_contour_filled(aes(y=c_l_contact_ratio,x=productivity_loss, z = economic_loss_ratio_ADJUSTED),size=0.009)+scale_fill_manual(values = wes_palette(14, name = "Zissou1", type = "continuous"), name = "Loss averted (%)")+theme_minimal()+ylab("Ratio of consumption-to-labor contacts") +xlab("Productivity loss (%)") +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_text(color="grey30"),legend.position = "none")+scale_y_continuous(breaks = seq(0,1.25,by = 0.25))+scale_x_continuous(breaks = seq(0,90,by = 10))
    #+geom_point(aes(x=14.45368, y=4.2877), fill="white", colour="#4f418b", shape=21,size = 1.5, stroke =1) 
    
    
#4c

 
 p1<-ggplot(Fig4a) +geom_raster(interpolate = T,aes(y=c_l_contact_ratio,x=productivity_loss, fill = round(cases_per_100k_ratio,3)))+scale_fill_gradientn(colours = wes_palette(14, name = "Zissou1", type = "continuous"), name = "Cases per 100k averted")+theme_minimal()+ylab("Ratio of consumption-to-labor contacts") +xlab("Productivity loss (%)") +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_text(color="grey30"),legend.position="none")
    
    
    
#4d     
 Fig4c<-read.csv("/~/Fig4/figure_4_ua_heatmap_data.csv", header=TRUE,sep=",") 
 
 p2<-ggplot(Fig4c) +
  geom_contour_filled(aes(y=unavoid_avoid_contact_ratio,x=productivity_loss, z = economic_loss_ratio_ADJUSTED),size=0.009)+scale_fill_manual(values = wes_palette(14, name = "Zissou1", type = "continuous"), name = "Loss averted (%)")+theme_minimal()+ylab("Ratio of unavoilable to avoidable contacts") +xlab("Productivity loss (%)") +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_text(color="grey30"),legend.position = "none")+scale_y_continuous(breaks = seq(0,1.25,by = 0.25))+scale_x_continuous(breaks = seq(0,90,by = 10)) 
  
  
  #4e
  
 
p3<-ggplot(Fig4c) +geom_raster(interpolate = T,aes(y=unavoid_avoid_contact_ratio,x=productivity_loss, fill = round(cases_per_100k_ratio,3)))+scale_fill_gradientn(colours = wes_palette(14, name = "Zissou1", type = "continuous"), name = "Cases per 100k averted")+theme_minimal()+ylab("Ratio of unavoilable to avoidable contacts") +xlab("Productivity loss (%)") +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_text(color="grey30"),legend.position="none") 
    
##new  with same breaks 
 ggplot(temperature[lev %in% c(1000, 300)], aes(lon, lat, z = air.z)) +
  geom_contour_fill(global.breaks = FALSE) +
  scale_fill_divergent() +
  facet_grid(~lev) 
  
 breaks<-read.csv("/Users/abento/Dropbox/Corona_epicon/FIGURES_MS/DATA/Fig4/NEW/4ac.csv", header=TRUE,sep=",")  
ggplot(breaks) +geom_contour_filled(aes(y=value,x=productivity_loss, z = economic_loss_ratio_ADJUSTED),size=0.009)+ facet_grid(~ratio) +scale_fill_manual(values = wes_palette(14, name = "Zissou1", type = "continuous"), name = "Loss averted (%)")+theme_minimal()+ylab("Ratio of unavoilable to avoidable contacts") +xlab("Productivity loss (%)") +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_text(color="grey30"))+scale_x_continuous(breaks = seq(0,90,by = 10)) +annotate("point", x = 14.45, y = 0.2798612, colour = "blue")  
  
  
    
    
#4f barplots
    
cases <- data.frame(pol=c("A", "B", "C"),
                value=c(1, 1,1))
p4<-ggplot(data=cases, aes(x=pol, y=value, fill=pol))+
  geom_bar(stat="identity")+scale_fill_manual(values=c("#BF812D","#C7EAE5","#01665E")) +scale_x_discrete(breaks=c("A","B","C"),
        labels=c("Convex","Linear","Concave"))+theme_minimal()+ylab("Total cases per 100,000 (ratio)")+xlab("Contact function shape") + theme_ipsum(axis_title_size = 24, base_size = 20)+theme(legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+coord_cartesian(ylim=c(0,1))+scale_y_continuous(breaks=seq(0,1, 0.25))
        
#4g barplots        
      
loss <- data.frame(pol=c("A", "B", "C"),
                value=c(0.95, 0.91,0.86))
p5<-ggplot(data=loss, aes(x=pol, y=value, fill=pol))+
  geom_bar(stat="identity")+scale_fill_manual(values=c("#BF812D","#C7EAE5","#01665E")) +scale_x_discrete(breaks=c("A","B","C"),
        labels=c("Convex","Linear","Concave"))+theme_minimal()+ylab("Individual loss averted (%)")+xlab("Contact function shape") + theme_ipsum(axis_title_size = 24, base_size = 20)+theme(legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank())+coord_cartesian(ylim=c(0,1))+scale_y_continuous(breaks=seq(0,1, 0.25))  
        
        
## assemblaged in illustrator




#######################################################################################################################        
#######################################################################################################################          
#######################################################################################################################          
#######################################################################################################################          
####################################         SI             ###########################################################  
#######################################################################################################################         

## lollipop for a combination

#libraries
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(gridExtra)
library(hrbrthemes)
test<-read.csv("/Users/GB/Dropbox/Corona_epicon/FIGURES_MS/DATA/SI/test.csv", header=TRUE,sep=",")

test1<-subset(test, measure=="cases per 100k ratio")


p<-ggplot(test1, aes(x=variable_value, y=measure_value)) + facet_grid(vars(variable_name), vars(measure),scales="free")+
   geom_segment( aes(x=variable_value, xend=variable_value, y=0.95, yend=measure_value), color="skyblue", alpha=0.7) +
   geom_point( fill="blue", size=2, colour="skyblue",pch=21) +
   theme_light() +
   coord_flip() +
   theme(
     panel.grid.major.y = element_blank(),
     panel.border = element_blank(),
     axis.ticks.y = element_blank()
   )+theme_minimal()+ylab("")+xlab("")+scale_y_continuous(breaks=seq(0.95,1,0.05))

test2<-subset(test, measure=="percentage of decentralized loss averted by targeting")          
p1<-ggplot(test2, aes(x=variable_value, y=measure_value)) + facet_grid(vars(variable_name), vars(measure),scales="free")+
  geom_segment( aes(x=variable_value, xend=variable_value, y=0.4, yend=measure_value), color="skyblue", alpha=0.7) +
  geom_point( fill="blue", size=2, colour="skyblue",pch=21) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )+theme_minimal()+ylab("")+xlab("")#+scale_y_continuous(breaks=seq(0.4,1.5,0.1)) 
  
  
  p<-ggplot(test, aes(x=variable_value, y=measure_value)) + facet_grid(vars(variable_name), vars(measure),scales="free")+
    geom_segment( aes(x=variable_value, xend=variable_value, y=0.4, yend=measure_value), color="skyblue", alpha=0.7) +
    geom_point( fill="blue", size=2, colour="skyblue",pch=21) +
   theme_bw() +
    coord_flip()
 p<-p+  theme(

      axis.ticks.y = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    )
 
 p<-p+theme(strip.background =element_rect(fill="white"))
 p<-p+ylab("")+xlab("")
 p
  
           
          
          