# Notes ------------------------------------------------------------------------

# This code takes in daily death and case data and seroprevalence data, and
# infers infection fatality rate for each age class

# Case and death data are reformatted in the "extend_dates.R" script and 
# seroprevalence data are reformatted in the "edit_serodata.R" script. The data
# folder already contains the reformatted results of these scripts, so they 
# do not need to be run again.

# Be sure to change this line to set the working directory to the folder 
# that you are working from 
setwd("~/R/sarscov2-ifr-nyc")

# Preamble ---------------------------------------------------------------------

# These are essential for the function of the code
library(tidyverse)
library(magrittr)
library(rstan)
library(foreach)
library(truncnorm)
options(mc.cores = 5)

# Clear environment
rm(list = ls())

# loads necessary functions
source("utils.R")

# Definitions and setup  -------------------------------------------------------

# suffix to differentiate between files
suffix <- "nyc"

# Death types: 'Confirmed' or 'Combined' 
# Combined includes both confirmed and probable
death_type <- "Confirmed"

# Number of draws
# The higher, the better and more accurate, but the longer it'll take to run
n_post <- 5000

# Redo the IFR fit?
redo_fit <- T

# Loading and formatting epidemiological data ----------------------------------

# NYC case and death age-stratified epidata: March 23 - May 17
# Extended by scaled daily case and death data to March 9 - May 17
# Age groups: 0-18, 18-44, 45-64, 65-74, 75+
# Source: NYC Dept. of Health 

# Why _extd? See "extend_dates.R". This code needs a sufficient number
# of age-specific epidata dates to run, and that script prepends 
# the available data by 14 days (derived from daily NYC case data).
epidata_source <- "./data/NYCHealth_cumul_data_extd.csv"

# Note: The getEpidata() function is very specific to parsing out the NYCHealth
# data. 
age_epidata <- getEpidata(death_type, epidata_source) 

# Organizing data
# Output contains columns with the following headings: 
# age_class, date, case_cumul, death_cumul
age_epidata <- age_epidata %>% 
  group_by(date, age_class, var) %>% 
  ungroup() %>% 
  arrange(age_class, date) %>% 
  pivot_wider(id_cols = c("age_class", "date"),
              names_from = "var",values_from = "Total") %>%
  mutate(date = as.Date(date))

# Creating seroprevalence draws ------------------------------------------------

# Source: Rosenberg et al., specifically Table 2, NYC, age groups
# "pop" column comes from the CUNY resource

# Why _edited? See "edit_serodata.R". The age classes from Rosenberg et al.
# have been adjusted to match up with the age classes from NYCHealth.
# 0-17 is set as "0" for now.
serovalue_data <- read_csv("./data/serodata_edited.csv")
pop_0to17 <- serovalue_data$pop[serovalue_data$age_class == "0-17"] #saving this value for later
  
serovalue_data <- serovalue_data %>%
  select(-"...1") %>%
  filter(age_class != "0-17") # because 0-17 draws get added later

# empty frame for bootstrapped draws
sero_draws <- data.frame(age_class = character(),
                         seropos = double(),
                         sim = integer(),
                         pop = double())

# creates draws for each age class using a truncated normal distribution
for (i in 1:length(serovalue_data$age_class)){
  sero_draws <- rbind(sero_draws, 
                      data.frame(age_class = serovalue_data$age_class[i], 
                                 seropos = rtruncnorm(n = n_post, a = 0, b=1,
                                 mean = serovalue_data$seropos_mu[i],
                                 sd = serovalue_data$seropos_sig[i]),
                                 sim = 1:n_post,
                                 pop = serovalue_data$pop[i]))
}
sero_draws[is.na(sero_draws)] = 0

# get seroprevalence for the 0-17 age class. This is done with the get0to17
# function, which estimates the 0-17 seroprevalence in NYC from 
# the age-specific seroprevalence trend observed in Spain, in a study done by
# Pollan et al. This study is the only one available to include fine age classes
# Source: Pollan et al.
trend_source <- "./data/serop_age_spain.csv"

sero_draws <- sero_draws %>%
  rbind(get0to17(trend_source = trend_source, n_post = n_post, 
         nyc_draws = sero_draws$seropos[sero_draws$age_class == "18-44"],
         pop = pop_0to17)) %>%
  arrange(age_class)

# Delay distributions ----------------------------------------------------------
# (One of) the next task: get rid of multiple serodatas if applicable


n_sero_weeks <- 1  # number of serosurvey weeks; NYC only has one,
                    # but this model could be used to assess more

# Filled in with mean and 95% CI from Rosenburg et al.
serodata <- tribble(
  ~date, ~median, ~q025, ~q975,
  "2020-04-23", 22.7, 21.5, 24.0
) %>% 
  mutate(date = as.Date(date))

# Previously chosen delays:
#   "report", 1.49, 0.065, NA, NA, 0.756, 0.0458, NA, NA,  # Scire et al. 2020
#   "report_death", 2.1, 0.055, NA, NA, 0.87,  0.039, NA, NA # DGS

# New delays
#  Sample collection to reporting: median 2, IQR 1-4, 90% 7
#  Diagnosis to death: median 8, IQR 4-16
delay_params <- tribble(
  ~delay, ~logmu, ~logmu.sd, ~logmu.low, ~logmu.high, ~logsigma,
  ~logsigma.sd, ~logsigma.low, ~logsigma.high,
  "inc", 1.57, NA, 1.44, 1.69, 0.65, NA, 0.56, 0.73,     # Bi et al. 2020
  "report", 0.693, 0.065, NA, NA, 1, 0.0458, NA, NA,  # Greene et al. 2021
  "symp_sero", 2.34, 0.114, NA, NA, 0.38, 0.26, NA, NA,  # Stringhini et al.2020
  "report_death", 2.08, 0.055, NA, NA, 1,  0.039, NA, NA # Thompson et al. 2020
) %>% 
  group_by(delay) %>% 
  mutate(
    logmu.sd = case_when(is.na(logmu.sd) ~ getSD(logmu, logmu.low, logmu.high),
                         T ~ logmu.sd),
    logsigma.sd = case_when(is.na(logsigma.sd)
                            ~ getSD(logsigma, logsigma.low, logsigma.high),
                            T ~ logsigma.sd)
  ) %>% 
  ungroup() %>% 
  mutate(mean = lnormMean(logmu, logsigma),   # Mean and SD on natural scale
         sd = lnormSD(logmu, logsigma))

nmax <- 100  # max number of days in PMFs

# Initialize PDF matrices
pdf_inc <- matrix(rep(0, nmax * n_post), nrow = n_post)
pdf_report <- pdf_inc
pdf_symp_sero <- pdf_inc
pdf_sero <- pdf_inc
pdf_report_death <- pdf_inc
pdf_symp_death <- pdf_inc
pdf_death <- pdf_inc

# Sample distributions
for (i in 1:n_post) {
  pdf_inc[i, ] <- setDelayPMF(delay_params, "inc", nmax, rnd = T)
  pdf_report[i, ] <- setDelayPMF(delay_params, "report", nmax, rnd = T)
  pdf_symp_sero[i, ] <- setDelayPMF(delay_params, "symp_sero", nmax, rnd = T)
  # convolve with incubation period:
  pdf_sero[i, ] <- convolve(pdf_inc[i, ],
                            rev(pdf_symp_sero[i, ]), type = "o")[1:nmax]
  pdf_report_death[i, ] <- setDelayPMF(delay_params, "report_death", nmax,
                                       rnd = T)
  # convolve with delay from symptoms to hospitalization:
  pdf_symp_death[i, ] <- convolve(pdf_report[i, ], rev(pdf_report_death[i, ]),
                                  type = "o")[1:nmax]
  # convolve with incubation period: 
  pdf_death[i, ] <- convolve(pdf_inc[i, ], rev(pdf_symp_death[i, ]),
                             type = "o")[1:nmax]  
}

# Inference --------------------------------------------------------------------
# Settings for stan
control <- list(adapt_delta = .9, max_treedepth = 12, metric = "dense_e")
nwarmup <- 4000
niter <- nwarmup + 1000
n_chains <- 5

# creates a vector to keep inferred infections
# Has a vector for each age class
inferred_daily_infections <- data.frame(infections = double())

# IFR Computation - repeats for every age class
IFRs <- tibble()
for (n in 1:length(unique(sero_draws$age_class))){
  
  agec = unique(sero_draws$age_class)[n] # current age class
  
  res_file <- paste0("age_stratified_IFRs_distsample_", agec, suffix, ".rds") 
  stanfit_file <- paste0("stanfit_age_stratified_IFRs_distsample_", agec,
                         suffix, ".rds") 
  
  if (!file.exists(res_file) | redo_fit) {
    
    ifr_stan <- stan_model("ifr_betabinomial.stan") # Compile stan model
    
    epidata <- filter(age_epidata, age_class == agec) # Cases and deaths per day
    
    pop_age <-  sero_draws$pop[sero_draws$age_class == agec][1] # Population of current age class
    
    ind_date <- which(epidata$date %in% serodata$date) # date index on which serosurvey was done

    deaths <- as.numeric(epidata$death_cumul[ind_date]) #  finds deaths on the serosurvey date

    n_data <- nrow(epidata) # finds how many days we have data

    
    # Thetas = serosurvey draws, in units of seroprevalent %, between 0 and 1
    # Formats posterior draws as matrix with each row corresponding to a 
    # serosurvey date, and each column to a single draw 
    thetas <- sero_draws %>% 
      filter(age_class == agec) %>% 
      select(seropos) %>% 
      t() %>% 
      as.matrix()
    
    # Isero = serosurvey applied to population, between 0 and the max population
    Isero <- pop_age * thetas
    
    # Creates empty matrices for Istar and phi [n data points x n samples]
    Istar <- matrix(rep(0, n_data * n_post), ncol = n_post)
    phi <- Istar
    
    # Loop over samples of delay distributions
    for (i in 1:n_post) {
      
      # Deconvolve state of cumulative infections up to proportionality (Istar)
      # Istar also equivalent to Isero times proportionality constant alpha
      # see Perez-Saez et al. supplemental material
      cases <- as.numeric(epidata$case_cumul)
      Istar[, i] <- computeIstar(cases, pdf_inc[i, ], pdf_report[i, ])
      
      # phi = fraction of infected population at risk of dying on day i
      phi[, i] <- (convolve(Istar[, i], rev(pdf_death[i, ]), type = "o")/
                     convolve(Istar[, i], rev(pdf_sero[i, ]),
                              type = "o"))[1:n_data]
    }
    
    # Selects the phi for date of serosurveys
    phis <- phi[ind_date, ]

    # I = number of Infected at risk of dying,  restricted to date of serosurvey
    I <- round(Isero * phis) 
    
    # Finds the daily infections from the cumulative infections
    # NOTE - the time series of inferred infections is given by Istar
    # Istar is in a cumulative form, and we need to get it to daily
    infections_agec <- vector()
    for (i in 2:n_data){
      I1 <- (mean(Istar[i-1,])) #infections for day before
      I2 <- (mean(Istar[i,])) #infections for day of
      
      infections_agec[i-1] <- I2 - I1 #finds daily infections
    }
    
    # binds it to the overall infections
    inferred_daily_infections <- rbind(inferred_daily_infections,
                                       infections_agec)
    
    # Set min of I to observed deaths to avoid errors in stan
    # Set max of I to total population in this age class 
    for (i in 1:nrow(I)) {
      I[i, I[i,] < deaths[i]] <- deaths[i] 
      I[i, I[i,] > pop_age] <- pop_age 
      I[i,] <- as.integer(I[i,])
    }
    
    # Data for stan model
    data <- list(N = n_sero_weeks, 
                 M = n_post, 
                 deaths = as.array(deaths), 
                 I = I)
    
    # Sample from stan model
    ifr_stanfit <- sampling(ifr_stan,
                            data = data,
                            init = rndInit, 
                            chains = n_chains,
                            warmup = nwarmup, 
                            iter = niter,
                            control = control,
                            save_warmup = FALSE)
    
    saveRDS(ifr_stanfit, file = "stanfit_file.rds")
    
    # Extract and save
    ifrs <- rstan::extract(ifr_stanfit, pars = "IFR")$IFR
    saveRDS(ifrs, file = "res_file.rds")
    
  } else {
    # Load results
    ifrs <- readRDS(res_file)
  }
  IFRs <- rbind(IFRs, tibble(ifr = ifrs, age_class = agec))
  }

# Export -----------------------------------------------------------------------

# Saves outputs to R Data Structure for postprocess
save(list = ls(), file = "data_setup.rda")
saveRDS(IFRs, file = paste0("age_stratified_IFRs_", suffix, ".rds"))
