# Edit serodata

# The goal of this script is to create an edited version of the serodata
# from Rosenberg et al. 2020, so that it's compatible with the main script.
#
# Edits include:
# 1) Manipulating the age classes so that they match those given by NYCHealth
#   Initial age classes: 18-44, 45-54, 55+
#   Ending age classes: 0-17, 18-44, 45-64, 65-74, 75+
# 2) Adding a column for population for each age class

# Clear environment
rm(list = ls())

library(tidyverse)
library(magrittr)

# Population by age 
# Contains the following columns with the following headings:
# age_class, age_midpoint, total
#
# Adapted from NYC population data, which is based on census data
# Source: 
#  https://www.baruch.cuny.edu/nycdata/population-geography/age_distribution.htm
age_popdata <- read.csv("./data/NYC_pop_age_total.csv")


# Source: Rosenberg et al., specifically Table 2, NYC, age groups
# "pop" column comes from the CUNY resource
# 18-34 and 35-44 were combined in a population-weighted average
serovalue_data <- read_csv("./data/serodata.csv") %>% 
  group_by(age_class)

# To get the 45-64 serodata mu and sig, some population weighting has to happen
serodata_45to54 <- c(serovalue_data[2,]$seropos_mu,
                     serovalue_data[2,]$seropos_sig)
serodata_55up <- c(serovalue_data[3,]$seropos_mu,
                   serovalue_data[3,]$seropos_sig)

# population-weight for the 45-54 serovalue = 
#               (45-54 population)/(45-64 population)
frac_45to54 <- serovalue_data[2,]$pop/
  age_popdata$total[age_popdata$age_class == "45-64"]

# In order: 0-17, 18-44, 45-64, 65-74, 75+
seropos_mu = c(0,
               serovalue_data$seropos_mu[1],
               frac_45to54*serovalue_data$seropos_mu[2] +
                 (1-frac_45to54)*serovalue_data$seropos_mu[3],
               serovalue_data$seropos_mu[3],
               serovalue_data$seropos_mu[3])
seropos_sig = c(0,
               serovalue_data$seropos_sig[1],
               frac_45to54*serovalue_data$seropos_sig[2] +
                 (1-frac_45to54)*serovalue_data$seropos_sig[3],
               serovalue_data$seropos_sig[3],
               serovalue_data$seropos_sig[3])

serodata <- data.frame(age_class = age_popdata$age_class,
                       seropos_mu = signif(seropos_mu,3),
                       seropos_sig = signif(seropos_sig,3),
                       pop = age_popdata$total)

write.csv(serodata, "./data/serodata_edited.csv")
