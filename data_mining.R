########################################################################
#                            COVID-19 DATA
########################################################################
# Last updated February 24, 2021 (Pacific Time)
# ***From: https://data.cdc.gov/Case-Surveillance/United-States-COVID-19-Cases-and-Deaths-by-State-o/9mfq-cb36
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm(list = ls())
pacman::p_load(tidyverse,lubridate)

# Directory
#setwd("/home/tatiane/Insync/tfr1@de.ufpe.br/Google Drive/RUBXII_regres_c19/GitHub_C19")
setwd("/media/tatiane/239fb286-e175-4578-b0d3-48c750946446/RUBXII_regres_c19/GitHub_C19")
data_CDC <- read_csv("United_States_COVID-19_Cases_and_Deaths_by_State_over_Time.csv")
#View(data_CDC)
sigla_state <- c(
  "AL", # - Alabama
  "AK", #  - Alaska
  "AZ", # - Arizona
  "AR", # - Arkansas
  "CA" , #- California
  "CO" , #- Colorado
  "CT" , #- Connecticut
  # "DC" , #- District of Columbia
  "DE" , #- Delaware
  "FL" , #- Florida
  "GA" , #- Georgia
  "HI" , #- Hawaii
  "ID" , #- Idaho
  "IL" , #- Illinois
  "IN" , #- Indiana
  "IA" , #- Iowa
  "KS" , #- Kansas
  "KY" , #- Kentucky
  "LA" , #- Louisiana
  "ME" , #- Maine
  "MD" , #- Maryland
  "MA" , #- Massachussets
  "MI" , #- Michigan
  "MN" , #- Minnesota
  "MS" , #- Mississippi
  "MO" , #- Missouri
  "MT" , #- Montana
  "NE" , #- Nebraska
  "NV" , #- Nevada
  "NH" , #- New Hampshire
  "NJ" , #- New Jersey
  "NM" , #- New Mexico
  "NY" , #- New York
  "NC" , #- North Carolina
  "ND" , #- North Dakota
  "OH" , #- Ohio
  "OK" , #- Oklahoma
  "OR", #- Oregon
  "PA" , #- Pennsylvania
  "RI" , #- Rhode Island
  "SC" , #- South Carolina
  "SD" , #- South Dakota
  "TN" , #- Tennessee
  "TX" , #- Texas
  "UT", # - Utah
  "VT", # - Vermont
  "VA", # - Virginia
  "WA", # - Washington
  "WV", # - West Virginia
  "WI", # - Wisconsin
  "WY" # - Wyoming
)

data_CDC1 <- data_CDC%>%dplyr::filter(state %in% sigla_state)%>%
  dplyr::mutate(as.factor(state))%>%
  dplyr::mutate(date = as.Date(format(as.Date(submission_date, '%m/%d/%Y'), '%Y/%m/%d')))%>%
  arrange(date)%>%
  dplyr::select(state,
                date,
                tot_cases,
                tot_death)

#View(data_CDC)

# First case's date  -->  Only for record 
data_CDC1_date <- data_CDC1 %>%
  dplyr::filter(tot_cases!=0) %>%
  dplyr::group_by(state) %>%
  dplyr::summarise(
    data_CDC1_date=min(date)
  )

#View(data_CDC1_date)  

# By selecting the total of deaths
panel <- data_CDC1 %>% 
  dplyr::filter(tot_cases>=10) %>%
  dplyr::group_by(state) %>%
  dplyr::slice(which.min(date)+30,
               which.min(date)+90,
               which.min(date)+180)%>%
  mutate(panel=1:3)  %>%
  dplyr::select(state, panel,
                date, deaths = tot_death) 

#View(panel)

#Covariates
cov <- read_csv("covariates.csv")

# By adding the covariates
panel1 <- panel%>%inner_join(cov,by="state")%>%
  mutate(MR = deaths/Pop_2020*100,
         T90 = ifelse(panel==2,1,0),
         T180 = ifelse(panel==3,1,0))%>%
  dplyr::select( State,
                 state,
                 panel,
                 MR,
                 PD,GINI,BEDS,SR,HDI,LE,PR,
                 T90,T180)%>%mutate(
                   panel = ifelse(panel==1,30,panel),
                   panel = ifelse(panel==2,90,panel),
                   panel = ifelse(panel==3,180,panel)
                 )%>%rename(TIME = panel)

#View(panel1)
#write_csv(panel1, "C19_MR_cov_USAbyState.csv")

