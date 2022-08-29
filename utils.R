# DESCRIPTION
# This file contains essential functions for the "draw_ifrs.R" file


#' @title Compute PDF or CDF of delay
#' @description Compute PDF or CDF of lognormal delays 
#' 
#' @param lmu mean on the logscale
#' @param lsigma sd on the logscale
#' 
#' @return vector of either probability densities or probabilities
#' 
computeDelay <- function(lmu, lsigma, nmax, dt = 1, start = dt/2, pdf = T) {
  if (pdf) {
    res <- dlnorm(seq(start, nmax, by = dt), meanlog = lmu, sdlog = lsigma)
  } else {
    res <- plnorm(seq(start, nmax, by = dt), meanlog = lmu, sdlog = lsigma)
  }
  return(res)
}

#' @title Get time delay PDF or CDF
#' @description Get the PDF or CDF of the delay distributions assuming a log-normal
#' 
#' @param delay_params Parameters of the log-normal distribution
#' @param delay Name of the delay
#' @param nmax Max number of timepoints, here using timestep dt
#' @param dt timestep in days along which to compute the densities
#' @param start firts timepoint to use
#' @param pdf boolean, whether to return the PDF or the CDF
#' @param rnd Whether to draw distribution parameter values
#' 
#' @return vector of either probability or probability density values
#' 
getDelay <- function(delay_params = NULL, delay = NULL, nmax, dt = 1, start = dt/2, pdf = T, rnd = F) {
  
  if (!(delay %in% delay_params$delay))
    stop("Delay name not known")
  
  ind <- which(delay_params$delay == delay) 
  
  if (!rnd) {
    # Set lognorm parameters mu and sigma on the logscale
    lmu <- delay_params$logmu[ind]
    lsigma <- delay_params$logsigma[ind]
  } else {
    # Set lognorm parameters mu and sigma on the logscale
    lmu <- truncnorm::rtruncnorm(1, a = 0, mean = delay_params$logmu[ind], sd = delay_params$logmu.sd[ind])
    lsigma <- truncnorm::rtruncnorm(1, a = 0, mean = delay_params$logsigma[ind], sd = delay_params$logsigma.sd[ind])
  }
  
  res <- computeDelay(lmu, lsigma, nmax, dt = dt, start = start, pdf = pdf)
  return(res)
}

#' @title Set Delay PMF
#' @description Set the PMF of the delay distributions assuming a log-normal
#' 
#' @param delay_params Parameters of the log-normal distribution
#' @param delay Name of the delay
#' @param nmax Max number of timepoints, here using a 1 day dt
#' @param rnd Whether to draw distribution parameter values
#' 
#' @details The daily PMF values are computed by taking the difference between the
#' delay distribution's CDF at times t and t-1: PMF(t) = CDF(t+1) - CDF(t).
#' 
#' @return vector of probability values
#' 
setDelayPMF <- function(delay_params, delay, nmax, dt = 1, rnd = F) {
  
  if (!(delay %in% delay_params$delay))
    stop("Delay name not known")
  
  ntot <- nmax/dt + 1 # total number of timesteps
  
  cdf <- getDelay(delay_params, delay, nmax, start = 0, pdf = F, rnd = rnd)
  
  prob_values <- cdf[2:ntot] - cdf[1:(ntot-1)]
  return(prob_values)
}

#' @title Compute Istar
#' @description Function to compute the deconvolution of reported cases and the 
#' incubation period accounting for reporting delay. Returns an estimate of the
#' cumulative number of infections up to a proportionality constant (Istar).
#'
#' @param cases time series of daily cumulative reported cases
#' @param pdf_inc PMF of the incubation period
#' @param pdf_report PMF of the reporting delay
#' @param gamma the threshold value for bounding the inverse fourier transform values
#'
#' @return A vector of the values of Istar
#' 
computeIstar <- function(cases, pdf_inc, pdf_report,  gamma = .05) {
  
  nc <- length(cases)
  ni <- length(pdf_inc)
  nr <- length(pdf_report)
  
  # FFT convolution from `convolve` function in R
  # Pad pdfs
  ntot <- ni + nr - 1
  pdf_inc2 <- c(pdf_inc, rep.int(0, ntot - ni))
  pdf_report2 <- c(pdf_report, rep.int(0, ntot - nr))
  ntot <- length(pdf_report2)
  F_pdf_comb <- fft(pdf_inc2) * fft(pdf_report2)
  pdf_comb <- Re(fft(F_pdf_comb, inverse  = T))/ntot
  # Preserve sum(prob) = 1
  pdf_comb <- pdf_comb/sum(pdf_comb)
  
  # Pad cases
  ntot2 <- nc + ntot - 1
  cases2 <- c(cases, rep.int(0, ntot2 - nc))
  eps <- 1e-10
  pdf_comb2 <- c(pdf_comb, rep(0, ntot2 - ntot))
  # fourier transform of convolution pdf
  F_cases <- fft(cases2)
  F_pdf_comb2 <- fft(pdf_comb2)
  
  # Water level regularization to prevent numerical instability
  # From https://cnx.org/resources/22c9f37591a06c51a0e632cc790ec83bcb853aa5/inverseFilter.m
  R <- F_pdf_comb2
  R1 <- F_pdf_comb2*(abs(F_pdf_comb2)>0)+1/gamma*(abs(F_pdf_comb2)==0)
  iR <- 1/F_pdf_comb2
  # invert the filter using threshold gamma
  G <- iR *((abs(R1)*gamma)>1)+gamma*abs(R1)*iR*((abs(R1)*gamma)<=1);
  Istar <- Re(fft(F_cases*G, inverse = T))[1:nc]
  # Sanity check for cumulative function
  Istar2 <- cummax(Istar)
  
  # # to get back to pdf
  # ntot <- length(pdf_comb2)
  # pdf_comb3 <- Re(fft(F_pdf_comb2, inverse  = T))/ntot
  # # Preserve sum(prob) = 1
  # pdf_comb3 <- pdf_comb3/sum(pdf_comb3)
  
  return(Istar2)
}

#' @title Compute Bayesian p-values
#' @description Compute the Bayesian p-value that two vectors of posterior draws 
#' have different means.
#' 
#' @param x vector of posterior draws of parameter to compare
#' @param y vector of posterior draws of reference parameter
#' @details the hypothesis that is tested is x = y using y - x = 0.
#' 
#' @return p-value
#' 
computePval <- function(x, y) {
  if (length(x) != length(y))
    stop("Vectors need to be the same length")
  
  return(min(sum((y-x) > 0), sum((y-x) < 0))*2/length(x))
}

#' @title Random initial values
#' @description Produces random initial values for Stan's HMC
#'
#' @return list with initial parameter values
#' 
rndInit <- function() {
  # Draw hyperprior parameters from priors 
  phi <- rbeta(1, 1, 6.5) # median of prior ~ 0.1
  lambda <- actuar::rpareto1(1, shape = 1.5, min = 0.1)
  # Draw IFR
  IFR <- runif(1, 1e-6, 1e-1)
  list(phi = phi, lambda = lambda, IFR = IFR)
}


#' @title Get SD
#' @description Computes the sd assuming a normal based on the 95% CI
#' 
#' @param mu mean
#' @param q025 0.025 quantile
#' @param q975 0.975 quantile
#' 
#' @return the standard deviation
#' 
getSD <- function(mu, q025, q975) {
  sd1 <- (mu - q025)/2
  sd2 <- (q975 - mu)/2
  return(mean(c(sd1, sd2)))
}

#' @title Log-norm mean
#' @description Computes the mean of a lognormal distribution
#' 
#' @param logmu mean on the logscale
#' @param logsigma sd on the logscale
#' 
#' @return the mean
#' 
lnormMean <- function(logmu, logsigma) {
  exp(logmu + .5 * logsigma^2)
}

#' @title Log-norm sd
#' @description Computes the sd of a lognormal distribution
#' @param logmu mean on the logscale
#' @param logsigma sd on the logscale
#' @return the sd
lnormSD <- function(logmu, logsigma) {
  sqrt((exp(logsigma^2) - 1) * exp(2 * logmu + logsigma^2))
}

#' @title Get quantiles
#' @description Computes the quantiles of a matrix of values
#' 
#' @param mat matrix over which to compute quantiles
#' 
#' @return dataframe of quantiles
#' 
getQuantiles <- function(mat) {
  mat %>% 
    as.data.frame() %>%
    mutate(sim = row_number()) %>%
    gather(var, value, -sim) %>%
    mutate(time = as.numeric(str_replace_all(var, "V", ""))) %>%
    group_by(sim) %>%
    arrange(time) %>%
    mutate(cumvalue = cumsum(value)) %>%
    group_by(time) %>%
    summarise(q025 = quantile(value, .025),
              q975 = quantile(value, .975),
              median = quantile(value, 0.5),
              mean = mean(value),
              cdf.q025 = quantile(cumvalue, .025),
              cdf.q975 = quantile(cumvalue, .975),
              cdf.median = quantile(cumvalue, 0.5),
              cdf.mean = mean(cumvalue))
}

#' @title organize epidata
#' @description subsets age-stratified cases and deaths from a dataset. Subsets
#' based on death type - Confirmed deaths or Confirmed+Probable Deaths. Renames
#' age classes based on whether we're adjusting the data
#' 
#' 
#' @return data frame of cumulative cases and deaths
#' 
getEpidata <- function(death_type, epidata_source) {

  if (death_type == "Confirmed"){
    epidata <- read.csv(epidata_source) %>%
      filter(var != "death_cumul_comb") %>%
      select(-X)
    epidata["var"][epidata["var"] == "death_cumul_conf"] <- "death_cumul"
  } else if (death_type == "Combined"){
    epidata <- read.csv(epidata_source) %>%
      filter(var != "death_cumul_conf") %>%
      select(-X)
    epidata["var"][epidata["var"] == "death_cumul_comb"] <- "death_cumul"
  }
  
  return(epidata)
}

#' @title get serosurvey draws for the 0-17 age class
#' @description Extrapolates 0-17 serosurvey draws from a source without
#' serodata for the 0-17 age class. Uses trend data from an external source
#' to generate a ratio between the 0-17 age class and another age class - here,
#' it's the 18-44 age class. Applies the ratio to the existing 18-44 data
#' 
#' @return data frame of serovalues for age 0-17
#' 
get0to17 <- function(trend_source, n_post, nyc_draws, pop){
  # SPECIAL CASE: Age class 0-17
  #   There is no 0-17 serodata available for New York City. However, there is
  #   0-17 case and death data available for New York City. We derived a scale of 
  #   0-17:18-45 seroprevalence from Spain data (Pollan et al. 2020) and then 
  #   applied that scale to the 18-44 seroprevalence in NYC to obtain the 0-17
  #   seroprevalence for NYC. We used Pollan et al. 2020 because it is the
  #   only sufficiently large study to include age-specific seroprevalence values
  #   in our target range
  
  # set up the ratio vector - each serosurvey draw gets its own scaling vector
  ratio_draws_0to17 <- rep(0, n_post)
  
  trend_data <- read.csv(trend_source)
  age_sero_fit = lm(mean~age_midpoint+I(age_midpoint^2), data = trend_data)
  trend_data$fitted = fitted(age_sero_fit)
  trend_data$scale = trend_data$fitted/mean(trend_data$mean)
  
  # The trend data is organized by the midpoints of each age class
  for (i in 1:n_post){
    
    # draws 1 value for each class between ages 0-17 of seroprevalence trend data
    nclasses_0to17 <- 
      length(trend_data$age_class[trend_data$age_midpoint < 18])
    trend_draws_0to17 <- rep(0, nclasses_0to17)
    for (j in 1:nclasses_0to17){
      trend_draws_0to17[j] <- rnorm(1,
                                    mean = (trend_data[j,]$mean)/100,
                                    sd = (trend_data[j,]$UpperCI-trend_data[j,]$mean)/100)
      if (trend_draws_0to17[j] < 0 | is.na(trend_draws_0to17[j])){
        trend_draws_0to17[j] <- 1e-6
      }
    }
    
    # draws 1 value for each class between ages 18-44 of seroprevalence trend data
    nclasses_18to44 <- 
      length(trend_data$age_class[trend_data$age_midpoint > 18 
                                  & trend_data$age_midpoint < 45])
    trend_draws_18to44 <- rep(0, nclasses_18to44)
    for (j in (nclasses_0to17+1):(nclasses_0to17+nclasses_18to44)){
      trend_draws_18to44[j-nclasses_0to17] <- rnorm(1, 
                                     mean = (trend_data[j,]$mean)/100, 
                                     sd = (trend_data[j,]$UpperCI-trend_data[j,]$mean)/100)
      if (trend_draws_18to44[j-nclasses_0to17] < 0 | is.na(trend_draws_18to44[j-nclasses_0to17])){
        trend_draws_18to44[j-nclasses_0to17] <- 1e-6
      }
    }
    
    # creates a single value for the 0-17:18-44 ratio out of the draws we 
    # made for this posterior iteration
    ratio_draws_0to17[i] <- mean(trend_draws_0to17)/mean(trend_draws_18to44)
  }
  print(mean(ratio_draws_0to17))
  
  # creates a 0-17 class from the ratio and the 18-44 seroprevalence draws from NYC
  draws_0to17 <- data.frame(age_class = "0-17", 
                            seropos = ratio_draws_0to17*nyc_draws,
                            sim = 1:n_post,
                            pop = pop)
  return(draws_0to17)
}
