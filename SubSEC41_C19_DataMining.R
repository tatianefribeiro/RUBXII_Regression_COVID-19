########################################################################
#                            COVID-19 DATA
########################################################################
# Last updated August 06, 2020 (Pacific Time)
# ***From: http://www.healthdata.org/covid/data-downloads
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm(list = ls())
library(tidyverse)

# Directory
setwd("Please, define here your directory")

ref_hosp <- read_csv("Reference_hospitalization_all_locs.csv")  # ***
ref_hosp <- ref_hosp%>%filter(location_id!=35)
cov <- read_csv("covariates.csv")

nstate = 50
name_state = c(
  "Alabama",
  "Alaska",
  "Arizona",
  "Arkansas",
  "California",
  "Colorado",
  "Connecticut",
  "Delaware",
  "Florida",
  "Georgia",
  "Hawaii",
  "Idaho",
  "Illinois",
  "Indiana",
  "Iowa",
  "Kansas",
  "Kentucky",
  "Louisiana",
  "Maine",
  "Maryland",
  "Massachusetts",
  "Michigan",
  "Minnesota",
  "Mississippi",
  "Missouri",
  "Montana",
  "Nebraska",
  "Nevada",
  "New Hampshire",
  "New Jersey",
  "New Mexico",
  "New York",
  "North Carolina",
  "North Dakota",
  "Ohio",
  "Oklahoma",
  "Oregon",
  "Pennsylvania",
  "Rhode Island",
  "South Carolina",
  "South Dakota",
  "Tennessee",
  "Texas",
  "Utah",
  "Vermont",
  "Virginia",
  "Washington",
  "West Virginia",
  "Wisconsin",
  "Wyoming"
)

ref_hosp1 <- ref_hosp%>%dplyr::filter(location_name %in% name_state)%>%
  dplyr::select(location_name,
                date,
                confirmed_infections,
                totdea_mean
  )%>%
  replace(is.na(.), 0)%>%
  dplyr::group_by(location_name) %>%
  mutate(confirmed_infections_cum = cumsum(confirmed_infections))  %>%
  ungroup()


# First case's date  -->  Only for record
case1_date <- ref_hosp1 %>% 
  dplyr::filter(confirmed_infections_cum!=0) %>%
  dplyr::group_by(location_name) %>%
  dplyr::summarize(
    case1_date=min(date)
  )

# By selecting the total of deaths
panel <- ref_hosp1 %>% 
  dplyr::filter(confirmed_infections_cum>=20) %>%
  dplyr::group_by(location_name) %>%
  slice(which.min(date)+30,
        which.min(date)+60,
        which.min(date)+90,
        which.min(date)+120)%>%
  mutate(panel=1:4)  %>%
  dplyr::select(location_name, panel,
                date, deaths = totdea_mean) 

which(panel[,4]==0)
panel[121,4] <- 44
panel[155,4] <- 833

## By adding the covariates
panel1 <- panel%>%inner_join(cov,by="location_name")%>%
  mutate(MR = deaths/Pop_2020*100,
         T60 = ifelse(panel==2,1,0),
         T90 = ifelse(panel==3,1,0),
         T120 = ifelse(panel==4,1,0))%>%
  dplyr::select( location_name,
                 panel,
                 MR,
                 PD, 
                 HDI,
                 GINI,
                 BEDS,
                 PR,
                 SR,
                 AT,
                 MA,
                 T60,
                 T90,
                 T120)

panel2 <- panel1%>%rename(STATE = location_name)%>%mutate(
  panel = ifelse(panel==1,30,panel),
  panel = ifelse(panel==2,60,panel),
  panel = ifelse(panel==3,90,panel),
  panel = ifelse(panel==4,120,panel)
)%>%rename(TIME = panel)

View(panel2)
write_csv(panel2, "COVID19_US_states.csv")

