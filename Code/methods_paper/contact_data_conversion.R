library(tidyverse)

setwd("../../Matrices_ForContact/USA/")

##### Numbers we use in the paper:
daily_cons_contacts = 4.2877 		# avg daily contacts at consumption activities
daily_work_contacts = 3.1625		# avg daily contacts at labor activities
daily_other_contacts = 5			# avg daily unavoidable contacts

contact_matrix <- matrix(c(daily_cons_contacts,0,0),c(0, daily_work_contacts, 0), c(0,0,daily_other_contacts))

##### HOW THE NUMBERS ABOVE WERE CALCULATED:
The spreadsheet here (https://www.dropbox.com/scl/fi/67e53n9jc4yapmaz8irl7/Aggregated.xlsx.xlsx?cloud_editor=gsheet&dl=0&rlkey=dznxdnxoqoz51p3npcwip0fop#gid=581298885) shows the calculation of age-weighted contact rates, based on taking colMeans of the contact matrices below (according to the main categorization, possibly). We weight the age-specific contacts by their prevalence in the population to get pop-weighted average numbers of contacts. We then compare that number of contacts against the sum(rowMeans(Mall)) (or something very like it), and back out the number of unavoidable "other" contacts as the residual between the total and daily_cons_contacts+daily_work_contacts.

##### Read in the contact matrix data
Mall<-read.csv("M.csv",sep=",", header=FALSE)
Mother<-read.csv("Mhome.csv",sep=",", header=FALSE)
Mschool<-read.csv("Mschool.csv",sep=",", header=FALSE)
Mwork<-read.csv("Mwork.csv",sep=",", header=FALSE)
Mconsumption<-read.csv("Mfun.csv",sep=",", header=FALSE)
Mlabor<-Mwork+Mschool # School has more students than teachers/staff, but we don't have separate matrices (and cross-contact matrices) for contacts between students and teachers/staff at schools. School is labor for teachers/staff, but consumption (really, investment) for students. We treat all school contacts as labor, rather than consumption. This is the more justifiable assumption, we think, because if teachers/staff don't show up students can't, while teachers/staff still go to schools when students aren't there.

## Sensitivity check -- how much do the numbers vary if we instead classify school as consumption?
Mlabor2 <- Mwork
Mconsumption2 <- Mconsumption + Mschool

gamma<-1/5.1 ##Infectious period

# Function to calculate the average number of contacts here, accounting for mixing, using the next-generation matrix approach
compute_avg_contacts <- function(R0,M){
    Eig = max(Re(eigen(M)$values))
    return(Eig)
}

## Contacts across all activities
### If we use the big matrices
compute_avg_contacts(2.6, Mall) # tau all 0.03288291  
sum(rowMeans(Mall)) 	# 
### If we use the sum of things from smaller matrices
compute_avg_contacts(2.6,Mother) + compute_avg_contacts(2.6,Mlabor) + compute_avg_contacts(2.6,Mconsumption)


# Contacts from sum of averages using next-generation method(?)
compute_avg_contacts(2.6,Mother) #tau other 0.1436656
compute_avg_contacts(2.6,Mlabor) #tau work 0.06785578
compute_avg_contacts(2.6,Mconsumption) #tau fun 0.09867632

sum(colMeans(Mother))
sum(colMeans(Mlabor))
sum(colMeans(Mconsumption))

compute_avg_contacts(2.6,Mlabor2) #tau work 0.06785578
compute_avg_contacts(2.6,Mconsumption2) #tau fun 0.09867632

# The implied tau from the economic structure (how many hours spent working, how many dollars consumed, etc) is back-calculated given the R0 from epicon_script.R

## computing R0
compute.R0 <- function(beta,M){
    Eig = max(Re(eigen(M)$values))
    return(Eig*beta/gamma)
}

## computing beta
compute_beta <- function(R0,M){
    Eig = max(Re(eigen(M)$values))
    return(R0*gamma/Eig)
}



