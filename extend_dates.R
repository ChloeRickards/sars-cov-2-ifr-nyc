# NOTE
# This script extends the cumulative data provided by NYCHealth from
# March 23-May 17 to March 9-Mary 17.
# Age-stratified cumulative case data is supplied for March 23-May 17 2020.
# Total (i.e. non-age-stratified) daily case data is supplied  March 1-June 15.

# In order to run the IFR program, we need data extending several days before
# March 23. We extrapolated age-stratified cumulative data by 14 days (March 6-
# March 22). To do this, we scaled the total daily case data by the 
# age-stratified fractions of cumulative cases for the first 14 of available
# data (March 23-April 6)

# SETUP -------------------------------------------------------------------
rm(list = ls())

daily_epidata_source <- read.csv("./data/NYCHealth_daily_data.csv")
cumul_epidata_source <- read.csv("./data/NYCHealth_cumul_data.csv")

age_classes <- unique(cumul_epidata_source$age_class)
n_age_classes <- length(age_classes)
ndays <- 14 # number of days to prepend


# Imports cumulative data, sets up  lists
# Source: NYC Department of Health
cases <- filter(cumul_epidata_source,var == "case_cumul")

conf_deaths <- filter(cumul_epidata_source, var == "death_cumul_conf")
comb_deaths <- filter(cumul_epidata_source, var == "death_cumul_comb")

# CASES------------------------------------------------------------------

# Gets total cases for each day
case_total <- vector()
for (i in 1:length(unique(cases$date))){
  case_total[i] <- sum(cases$Total[cases$date == unique(cases$date)[i]])
}

# Gets fractions of cases for first 14 days
case_fractions <- vector()
for (i in 1:ndays){
  for(j in n_age_classes:1){
    case_fractions[length(case_fractions) + 1] <- 
      cases$Total[n_age_classes*i-(j-1)]/case_total[i]
  }
}
case_fractions <- data.frame(case_fractions, age_classes)

# Creates the scale we'll apply to the daily cases
case_scale <- data.frame(frac=rep(0,n_age_classes), age_classes)
for (i in 1:n_age_classes){
  frac_class <- filter(case_fractions, age_classes == age_classes[i])
  case_scale$frac[i] <- sum(frac_class$case_fractions)/ndays
}

# Retrieves daily case data from before March 23 and then scales the case data
# scales the daily case data down so that the cumulative cases from 
# March 6-22 matches the total number of cumulative cases on March 23
# (i.e. that the two sources become placed on the same scale)
# Note: the initial discrepancy is likely due to the difference between
# confirmed and probable cases.
# Note 2: The last line of this chunk removes the last row of dailies so that
# the dailies ends on March 22. We do this because the age-stratified
# cumulative cases start on March 23.

day1 <- as.Date(cumul_epidata_source[1,]$date)

cdailies <- daily_epidata_source %>%
  filter(var == "case_daily") %>%
  filter(date < day1+1) %>%
  filter(date > day1-ndays-1)
cdailies$Total <- cdailies$Total*(case_total[1]/sum(cdailies$Total))
cdailies <- head(cdailies,-1)

# This converts the total daily cases into cumulative age-stratified daily cases
# Epidata columns: age_class, Total, var, date
case_epidata <- data.frame()
for (i in 1:length(cdailies[,1])){
  if (i == 1){
    cumul_cases <- round(cdailies[i,]$Total*case_scale$frac)
  } else{
    cumul_cases <- prev_cases + round(cdailies[i,]$Total*case_scale$frac)
  }
  prev_cases <- cumul_cases
  
  for (j in 1:n_age_classes){
    case_epidata <- rbind(case_epidata, c(age_class = age_classes[j],
                                          Total = cumul_cases[j],
                                          var = "case_cumul",
                                          date = cdailies$date[i]))
  }
}



# DEATHS --------------------------------------------------------------
# Same-ish process as above, just repeated for the deaths
# Some variations due to how the death data is organized

# Gets total combined deaths
comb_death_total <- vector()
for (i in 1:length(unique(comb_deaths$date))){
  comb_death_total[i] <- 
    sum(comb_deaths$Total[comb_deaths$date == unique(comb_deaths$date)[i]])
}

conf_death_total <- vector()
for (i in 1:length(unique(conf_deaths$date))){
  conf_death_total[i] <- 
    sum(conf_deaths$Total[conf_deaths$date == unique(conf_deaths$date)[i]])
}

# Gets fractions of cases for first 14 days
comb_death_fractions <- vector()
for (i in 1:ndays){
  for(j in n_age_classes:1){
    comb_death_fractions[length(comb_death_fractions) + 1] <- 
      comb_deaths$Total[n_age_classes*i-(j-1)]/comb_death_total[i]
  }
}
comb_death_fractions <- data.frame(comb_death_fractions, age_classes)

# Creates the scale we'll apply to the daily cases
comb_death_scale <- data.frame(frac=rep(0,n_age_classes), age_classes)
for (i in 1:n_age_classes){
  frac_class <- filter(comb_death_fractions, age_classes == age_classes[i])
  comb_death_scale$frac[i] <- sum(frac_class$comb_death_fractions)/ndays
}

# Gets fractions of cases for first 14 days
conf_death_fractions <- vector()
for (i in 1:ndays){
  for(j in n_age_classes:1){
    conf_death_fractions[length(conf_death_fractions) + 1] <- 
      conf_deaths$Total[n_age_classes*i-(j-1)]/conf_death_total[i]
  }
}
conf_death_fractions <- data.frame(conf_death_fractions, age_classes)

# Creates the scale we'll apply to the daily cases
conf_death_scale <- data.frame(frac=rep(0,n_age_classes), age_classes)
for (i in 1:n_age_classes){
  frac_class <- filter(conf_death_fractions, age_classes == age_classes[i])
  conf_death_scale$frac[i] <- sum(frac_class$conf_death_fractions)/ndays
}

# Loads daily death data
comb_ddailies <- daily_epidata_source %>%
  filter(var == "death_daily_comb") %>%
  filter(date < day1+1) %>%
  filter(date > day1-ndays-1)
comb_ddailies$Total <- comb_ddailies$Total*(comb_death_total[1]/
                                            sum(comb_ddailies$Total))
comb_ddailies <- head(comb_ddailies,-1)

conf_ddailies <- daily_epidata_source %>%
  filter(var == "death_daily_conf") %>%
  filter(date < day1+1) %>%
  filter(date > day1-ndays-1)
conf_ddailies$Total <- conf_ddailies$Total*(conf_death_total[1]/
                                            sum(conf_ddailies$Total))
conf_ddailies <- head(conf_ddailies,-1)

# Converts the total daily deaths into cumulative age-stratified daily deaths
# Epidata columns: age_class, Total, var, date
comb_death_epidata <- data.frame()
for (i in 1:length(comb_ddailies[,1])){
  if (i == 1){
    cumul_deaths <- round(comb_ddailies[i,]$Total*comb_death_scale$frac)
  } else{
    cumul_deaths <- prev_deaths + 
      round(comb_ddailies[i,]$Total*comb_death_scale$frac)
  }
  prev_deaths <- cumul_deaths
  
  for (j in 1:n_age_classes){
    comb_death_epidata <- rbind(comb_death_epidata, c(age_class = age_classes[j],
                                            Total = cumul_deaths[j],
                                            var = "death_cumul_comb",
                                            date = comb_ddailies$date[i]))
  }
}

conf_death_epidata <- data.frame()
for (i in 1:length(conf_ddailies[,1])){
  if (i == 1){
    cumul_deaths <- round(conf_ddailies[i,]$Total*conf_death_scale$frac)
  } else{
    cumul_deaths <- prev_deaths + 
      round(conf_ddailies[i,]$Total*conf_death_scale$frac)
  }
  prev_deaths <- cumul_deaths
  
  for (j in 1:n_age_classes){
    conf_death_epidata <- rbind(conf_death_epidata, c(age_class = age_classes[j],
                                                      Total = cumul_deaths[j],
                                                      var = "death_cumul_conf",
                                                      date = conf_ddailies$date[i]))
  }
}

# OUTPUT ------------------------------------------------------

# Updating column names for the two dataframes we made
colnames(case_epidata) <- c("age_class", "Total", "var", "date")
colnames(comb_death_epidata) <- c("age_class", "Total", "var", "date")
colnames(conf_death_epidata) <- c("age_class", "Total", "var", "date")

# Duplicating the epidata so that we get a 

epidata <- rbind(conf_death_epidata, conf_deaths, comb_death_epidata,
                 comb_deaths, case_epidata, cases)

# If we're adjusting: duplicating the epidata so that we get values for both
# raw and adjusted age classes.

write.csv(epidata, "./data/NYCHealth_cumul_data_extd.csv")

