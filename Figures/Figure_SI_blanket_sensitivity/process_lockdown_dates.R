# Script to process covid lockdown dates and get mean duration and mean start date after first case in the US

library(tidyverse)
library(lubridate)

data <- read_csv("COVIDpolicies_2020-08-16.csv")

first_case <- "20-Jan-20"

df <- data %>% 
      mutate(interval_case1_to_school_close = interval(dmy(first_case),dmy(DateSchoolClose) )) %>%
      mutate(duration_till_school_close = time_length(interval_case1_to_school_close, unit="days")) %>%
      mutate(interval_case1_to_restaurant_close = interval(dmy(first_case),dmy(DateRestaurantOtherRestrict) )) %>%
      mutate(duration_till_restaurant_close = time_length(interval_case1_to_restaurant_close, unit="days")) %>%
      mutate(duration_till_first_close = pmin(duration_till_school_close,duration_till_restaurant_close,na.rm=TRUE)) %>%
      mutate(full_duration_till_first_reopen = time_length(interval(dmy(first_case),dmy(DateReOpen)),unit="days")) %>%
      mutate(full_duration_till_second_reopen = time_length(interval(dmy(first_case),dmy(DateReOpenPhase2)),unit="days")) %>%
      mutate(duration_till_first_reopen = full_duration_till_first_reopen - duration_till_first_close) %>%
      mutate(duration_till_second_reopen = full_duration_till_second_reopen - duration_till_first_close)


df$duration_till_first_close
mean(df$duration_till_first_close,na.rm=TRUE)
df$duration_till_first_reopen
mean(df$duration_till_first_reopen,na.rm=TRUE)

write.csv(df,file="lockdowns.csv")
