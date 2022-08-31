# Notes ------------------------------------------------------------------------
# This code requires the "data_setup.rda" file, and 
# "age_stratified_IFRs_[suffix].rds" which are created from running
# "draw_ifrs.R". The default rda and rds files contain results from
# the NYC analysis using confirmed deaths only.
#
# This produces all of the figures and tables from the paper

# Preamble and Setup -----------------------------------------------------------
library(tidyverse)
library(magrittr)
library(gridExtra)
library(zoo)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(scales)
library(betareg)
library(effects)

source("utils.R")

# Loading outputs from draw_ifrs
load("data_setup.rda")


# IFR Table --------------------------------------------------------------------
# Creates table with Age Class, Population, Deaths, Infection Prevalence,
# Estimated Number Infected, IFRs 


# Load and format IFRs
IFRs <- readRDS(paste0("age_stratified_IFRs_", suffix, ".rds"))
age_classes <- unique(sero_draws$age_class)
IFRs <- IFRs %>%
  group_by(age_class) %>%
  mutate(sim = row_number()) %>%
  ungroup() %>%
  inner_join(sero_draws,
             by = c("sim", "age_class" = "age_class"))

# Add seroprevelance and population data
IFR_stats <- 
  # Epi data
  age_epidata %>%
  group_by(age_class) %>% 
  arrange(desc(case_cumul)) %>% 
  select(case_cumul, death_cumul) %>% 
  slice(1) %>% ungroup() %>% 
  cbind(total=unique(sero_draws$pop)) %>%
  # Seroprevalence estimates
  inner_join(sero_draws %>% 
               group_by(age_class) %>% 
               summarise(sero.mean = mean(seropos),
                         sero.025 = quantile(seropos, 0.025),
                         sero.975 = quantile(seropos, 0.975))) %>% 
  inner_join(sero_draws %>% 
               group_by(age_class) %>% 
               summarise(infected.mean = mean(seropos),
                         infected.025 = quantile(seropos, 0.025),
                         infected.975 = quantile(seropos, 0.975))) %>%
  mutate(infected.mean = infected.mean*total,
         infected.025 = infected.025*total,
         infected.975 = infected.975*total) %>%
  # IFR estimates
  inner_join(
    IFRs %>% 
      group_by(age_class) %>% 
      mutate(ifr = ifr*1e2) %>% 
      summarise(mean = mean(ifr),
                q025 = quantile(ifr, 0.025),
                q975 = quantile(ifr, 0.975))
  ) %>% 
  group_by(age_class) %>% 
  mutate_at(vars(mean, q025, q975), 
            function(x) format(x, digits = 2)) %>% 
  ungroup() %>% 
  mutate(ifr = paste0(mean, " (", q025, "-", q975, ")"),
         sero = paste0(format(signif(100*sero.mean,3), big.mark = ","),
                       " (", format(signif(100*sero.025,3), big.mark = ","), 
                       "-", 
                       format(signif(100*sero.975,3), big.mark = ","), ")"),
         infected = paste0(format(signif(infected.mean,3), big.mark = ","),
                           " (", format(signif(infected.025,3), big.mark = ","), 
                           "-", 
                           format(signif(infected.975,3), big.mark = ","), ")"),
         age_class = factor(age_class, 
                            levels = age_classes[c(1, 2, 3,
                                                   4:length(age_classes))])) %>%
  arrange(age_class)

# Create and format the results table
IFR_table <- IFR_stats %>%
  select(age_class, total, death_cumul, sero, infected, ifr) %>% 
  rename(`Age class` = age_class,
         `Population` = total,
         `Deaths` = death_cumul,
         `Estimated Infection Prevalence` = sero,
         `Estimated Infected` = infected,
         IFR = ifr)

write.csv(IFR_table, file = paste0("./results/IFRs_table_", suffix, ".csv"))

# Case, death, and infection timeline ------------------------------------------

# Based on confirmed daily cases
# Dates of interest: 2020-03-01 to 2020-06-15
# Source: NYC Dept. of Health and Mental Hygiene 
# https://www1.nyc.gov/site/doh/covid/covid-19-data-archive.page 

# Formatting daily case data
daily_cases <- read_csv("./data/NYCHealth_daily_data.csv") %>%
  filter(var == "case_daily") %>%
  mutate(var = "Daily cases (7-day average)") %>%
  rename(DailyCount = Total)

avg_cases = rollmean(daily_cases$DailyCount, 7)
daily_cases <- daily_cases %>%
  tail(length(avg_cases)) %>%
  mutate(DailyCount = avg_cases)

# Formatting daily death data
if(death_type == "Confirmed"){
  daily_deaths <- read_csv("./data/NYCHealth_daily_data.csv") %>%
    filter(var == "death_daily_conf") %>%
    mutate(var = "Daily deaths") %>%
    rename(DailyCount = Total)
} else{
  daily_deaths <- read_csv("./data/NYCHealth_daily_data.csv") %>%
    filter(var == "death_daily_comb") %>%
    mutate(var = "Daily deaths") %>%
    rename(DailyCount = Total)
}

# Serosurvey midpoint, marked as a vertical line
serosurvey_midpoint <- as.Date("2020-04-23")

# Formatting infections
# Pulls infections from the Istar calculated in previous script
inferred_daily_infections <- distinct(inferred_daily_infections)
daily_infections <- vector()
for (i in 1:length(inferred_daily_infections)){
  daily_infections[i] <- sum(inferred_daily_infections[,i])
}
daily_infections <- rollmean(daily_infections[daily_infections > 0],7)/10
daily_infections <- data.frame(DailyCount = daily_infections,
                               date = head(epidata$date,
                                           length(daily_infections)),
                               var = "Inferred daily infections/10 (7-day average)")
daily_infections <- daily_infections[daily_infections$date <= serosurvey_midpoint,]


# Combines all of the relevant dates
timeline <- daily_deaths %>%
  rbind(daily_cases) %>%
  rbind(daily_infections) %>%
  rbind(data.frame(DailyCount=c(0,7500),
                   date=serosurvey_midpoint,
                   var = "Serosurvey midpoint"))%>%
  rename(Legend = var)

# Plot
ggplot(timeline, aes(x=date,y = DailyCount)) +
  geom_rect(aes(xmin=as.Date("2020-04-19"),
                xmax=as.Date("2020-04-28"),
                ymin=0,
                ymax=7500,
                fill="Serosurvey"),
            alpha=0.5)+
  scale_fill_manual('',
                    values = 'grey',
                    guide = guide_legend(override.aes = list(alpha = 1)))+
  #scale_color_manual(values=line$color)+
  geom_line(aes(color = Legend, size = Legend), size = 1.5)+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  theme_bw()+
  theme(axis.title=element_text(size=25),legend.position = c(.77, .8),
        axis.text=element_text(size=20),legend.text=element_text(size=18),
        legend.title=element_blank())+
  scale_size_manual(values = c(1, 1, 1, 1, 3))+
  scale_color_manual(values=c("#61ddff", "#2051e3", "#de2509", "black"))+
  scale_y_continuous(breaks = seq(0,8000,1000),labels=c("0","1,000","2,000",
                                                        "3,000","4,000","5,000",
                                                        "6,000","7,000","8,000"))+
  xlab("Date")+
  ylab("COVID-19 cases and deaths\n")

ggsave("./results/nyc_timeline.png", width = 13.5, height = 10.3)


# IFR Comparison ---------------------------------------------------------------

# Loads the IFRs to compare
# Organized by age midpoint, author, Region, mean IFR, lower CI, upper CI
IFR_comparison <- read_csv("./data/IFR_comparison.csv") %>%
  rename(Region = "DatasetCountry",
         LowerCI = "Lower CI",
         UpperCI = "Upper CI")


# Midpoints for our age classes: 0-18, 18-45, 45-65, 65-75, 75+
# Formatting our data to go with the comparison figure
NYC_compare <- data.frame(AgeMidpoint = c(9, 32, 55, 70, 85),
                          DatasetAuthor = "This study",
                          Region = "New York City",
                          IFR = c(IFR_stats$mean[IFR_stats$age_class == "0-17"],
                                  IFR_stats$mean[IFR_stats$age_class == "18-44"],
                                  IFR_stats$mean[IFR_stats$age_class == "45-64"],
                                  IFR_stats$mean[IFR_stats$age_class == "65-74"],
                                  IFR_stats$mean[IFR_stats$age_class == "75+"]),
                          LowerCI = c(IFR_stats$q025[IFR_stats$age_class == "0-17"],
                                      IFR_stats$q025[IFR_stats$age_class == "18-44"],
                                      IFR_stats$q025[IFR_stats$age_class == "45-64"],
                                      IFR_stats$q025[IFR_stats$age_class == "65-74"],
                                      IFR_stats$q025[IFR_stats$age_class == "75+"]),
                          UpperCI = c(IFR_stats$q975[IFR_stats$age_class == "0-17"],
                                      IFR_stats$q975[IFR_stats$age_class == "18-44"],
                                      IFR_stats$q975[IFR_stats$age_class == "45-64"],
                                      IFR_stats$q975[IFR_stats$age_class == "65-74"],
                                      IFR_stats$q975[IFR_stats$age_class == "75+"]))

# Stitching all the IFRs together to compare
IFR_comparison <- rbind(IFR_comparison, NYC_compare) %>%
  mutate(IFR = as.numeric(IFR),
         LowerCI = as.numeric(LowerCI),
         UpperCI = as.numeric(UpperCI))

# Formatting
options(scipen=10000)

# Plot on log scale
ggplot(IFR_comparison, aes(y=IFR, x=AgeMidpoint,col=Region)) +
  geom_line(aes(group = Region),position=position_dodge(width=2))+
  theme_bw()+
  geom_point(aes(shape=Region),size=3.5,position=position_dodge(width=2))+
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI),
                width = 3,position=position_dodge(width=2))+
  scale_y_log10(breaks = c(0.0001,0.001,0.01,0.1,1,10),labels=c("0%","0.001%","0.01%","0.1%","1.0%","10.0%"))+
  scale_x_continuous(limits=c(0,100),breaks=c(20,40,60,80,100))+
  theme(axis.title=element_text(size=25),legend.position = c(.2, .80),
        axis.text=element_text(size=25),legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.margin = margin(10,30,10,10))+
  scale_shape_manual(values=c(15, 16, 17,
                              15, 16, 17,
                              15, 16))+
  scale_color_manual(values=c("#abbf13",  "#24b379", "#000000",
                              "#d97b09", "#cc0000","#0c09d9",
                              "#ad05f5", "#f542e3"))+
  xlab("Age")+
  ylab("Infection fatality rate (log scale)")

ggsave("./results/ifr_comparison.png", width = 10, height = 10)

# Plot on linear scale
IFR_comparison$UpperCI[IFR_comparison$Region == "Sweden" & IFR_comparison$AgeMidpoint > 88] = 20
ggplot(IFR_comparison, aes(y=IFR, x=AgeMidpoint,col=Region)) +
  geom_line(aes(group = Region),position=position_dodge(width=2))+
  theme_bw()+
  geom_point(aes(shape=Region),size=4,position=position_dodge(width=2))+
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI),
                width = 0.25,position=position_dodge(width=2))+
  scale_y_continuous(breaks = c(0,5,10,15,20),labels=c("0%","5%","10%","15%","20%"), limits=c(0,20))+
  theme(axis.title=element_text(size=25),legend.position = c(.2, .80),
        axis.text=element_text(size=25),legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        plot.margin = margin(10,30,10,10))+
  scale_shape_manual(values=c(15, 16, 17,
                              15, 16, 17,
                              15, 16))+
  scale_color_manual(values=c("#abbf13",  "#24b379", "#000000",
                              "#d97b09", "#cc0000","#0c09d9",
                              "#ad05f5", "#f542e3"))+
  xlab("Age")+
  ylab("Infection fatality rate")

ggsave("./results/ifr_comparison_untransformed.png", width = 10, height = 10)



# Delay distribution plot ------------------------------------------------------

# Delays for plotting
nmax <- 50  # max number of days in PDFS
dt <- .1
ntot <- nmax/dt
n_dist_samples <- 1000

# Initialize PDF matrices
pdf_inc <- matrix(rep(0, ntot * n_dist_samples), nrow = n_dist_samples)
pdf_report <- pdf_inc
pdf_case <- pdf_inc
pdf_symp_sero <- pdf_inc
pdf_sero <- pdf_inc
pdf_report_death <- pdf_inc
pdf_symp_death <- pdf_inc
pdf_death <- pdf_inc

# Sample distributions
for (i in 1:n_dist_samples) {
  pdf_inc[i, ] <- getDelay(delay_params, "inc", nmax, dt = dt, rnd = T)*dt
  pdf_report[i, ] <- getDelay(delay_params, "report", nmax, dt = dt, rnd = T)*dt
  pdf_case[i, ] <- convolve( pdf_inc[i, ], rev(pdf_report[i, ]), type = "o")[1:ntot]
  pdf_symp_sero[i, ] <- getDelay(delay_params, "symp_sero", nmax, dt = dt, rnd = T)*dt
  
  # convolve with incubation period
  pdf_sero[i, ] <- convolve(pdf_inc[i, ], rev(pdf_symp_sero[i, ]), type = "o")[1:ntot]  
  pdf_report_death[i, ] <- getDelay(delay_params, "report_death", nmax, dt = dt, rnd = T)*dt
  
  # convolve with delay from symptoms to hospitalization
  pdf_symp_death[i, ] <- convolve(pdf_report[i, ],
                                  rev(pdf_report_death[i, ]),
                                  type = "o")[1:ntot] 
  
  # convolve with incubation period
  pdf_death[i, ] <- convolve(pdf_inc[i, ], rev(pdf_symp_death[i, ]), type = "o")[1:ntot]
}

# Initializing labels
delay_dict <- c(
  "i2s" = "Incubation period \n (Infection to symptom onset)",
  "s2c" = "Symptom onset to case reporting",
  "s2sero" = "Symptom onset to seroconversion",
  "c2d" = "Case reporting to death",
  "i2c" = "Infection to case reporting",
  "i2sero" = "Infection to seroconversion",
  "i2d" = "Infection to death"
)

legend_title <- "Delay distribution"

# Put together pdfs and cdfs
delays_2plot <- map2_dfr(
  list(pdf_inc, pdf_report, pdf_symp_sero, pdf_report_death, pdf_case, pdf_sero,
       pdf_death),
  c("i2s", "s2c", "s2sero", "c2d", "i2c", "i2sero", "i2d"), 
  ~getQuantiles(.x) %>%
    mutate(days = seq(dt/2, nmax, by = dt),
           var = .y)) %>% 
  mutate(type = case_when(
    var %in% c("i2s", "s2c", "s2sero", "c2d") ~ "literature",
    T ~ "derived"),
    var_name = delay_dict[var]) 

# CDFs for literature-derived delays
p_lit <- delays_2plot %>% 
  filter(type == "literature") %>% 
  ggplot(aes(x = days, y = cdf.median, ymin = cdf.q025, ymax = cdf.q975,
             fill = var_name, lty = var_name)) + 
  geom_ribbon(alpha = .25) +
  geom_line(aes(color = var_name)) +
  guides(color = guide_legend(legend_title),
         fill = guide_legend(legend_title),
         lty = guide_legend(legend_title)) +
  theme_bw() +
  labs(x = "time [days]", y = "CDF") +
  ggthemes::scale_color_few() +
  ggthemes::scale_fill_few() +
  theme(legend.position = c(.71, .2))

# CDFs for convolution-derived delays
p_conv <- delays_2plot %>% 
  filter(type == "derived") %>% 
  ggplot(aes(x = days, y = cdf.median, ymin = cdf.q025, ymax = cdf.q975,
             fill = var_name, lty = var_name)) + 
  geom_ribbon(alpha = .25) +
  geom_line(aes(color = var_name)) +
  guides(color = guide_legend(legend_title),
         fill = guide_legend(legend_title),
         lty = guide_legend(legend_title)) +
  theme_bw() +
  labs(x = "time [days]", y = "CDF") +
  theme(legend.position = c(.71, .2))

# Plot
p_comb <- arrangeGrob(grobs = list(p_lit, p_conv), nrow = 1)
ggsave(plot = p_comb, filename = "./results/delay_dist.png", width = 9, height = 4.5)


# Plotting PDFs
pdfs <- delays_2plot %>% 
  filter(var == "i2s" | var == "i2c" | var== "i2sero" | var== "i2d") %>% 
  ggplot(aes(x = days, y = median, ymin = q025, ymax = q975,
             fill = var_name, lty = var_name)) + 
  geom_ribbon(alpha = .25) +
  geom_line(aes(color = var_name)) +
  guides(color = guide_legend(legend_title),
         fill = guide_legend(legend_title),
         lty = guide_legend(legend_title)) +
  theme_bw() +
  labs(x = "time [days]", y = "PDF") +
  scale_color_manual(values = c("orange", "seagreen", "red","blue"))+
  scale_fill_manual(values= c("coral", "seagreen1", "firebrick", "lightblue"))+
  theme(legend.position = c(.71, .7),
        axis.title=element_text(size=20),
        axis.text=element_text(size=15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))

p_comb <- arrangeGrob(grobs = list(pdfs), nrow = 1)
plot(p_comb)



# Delay distribution violin plot -----------------------------------------------

# Delays of interest: here we are measuring infection delays
delay_list <- c("i2s", #infection to symptom
                "i2c",  #infection to case report
                "i2sero", #infection to seroconversion
                "i2d")  #infection to death
delay_names <- c("Incubation period \n (Infection to symptom onset)",
                 "Infection to case reporting",
                 "Infection to seroconversion",
                 "Infection to death")

# Going from PDFs to days distribution
n <- 10000
violin_data <- data.frame()
for (i in 1:length(unique(delays_2plot$days))){
  for (j in 1:length(delay_list)){
    if (round(delays_2plot$mean[delays_2plot$var == delay_list[j]][i]*n) > 0){
      violin_data <- rbind(violin_data,
                           data.frame(days = rep(delays_2plot$days[i],
                                                 times = 
                                                   round(delays_2plot$mean[delays_2plot$var==delay_list[j]][i]*n)),
                                      delay = delay_names[j]))
    }
  }
}

# Establish the box and whisker points, based on the median, IQR, and 95% CI
quantiles_95 <- function(x) {
  box_limits <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(box_limits) <- c("ymin", "lower", "middle", "upper", "ymax")
  box_limits
}

# Plot
ggplot(violin_data, aes(x=delay, y=days, fill=delay)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("#cc000080", "#d97b0980", "#24b37980","#0c09d980"))+
  stat_summary(fun.data = quantiles_95, geom="boxplot", fill = "white", width = 0.1, lwd = 1) +
  scale_y_continuous(breaks=seq(0,50,5))+
  theme_bw()+
  theme(axis.title=element_text(size=20), axis.text=element_text(size=17),
        legend.position = "none")+
  xlab("Delay distribution")+
  ylab("Days")

ggsave("./results/delay_violins.png", width = 12, height = 8)


# Gini and Income vs. IFR plots ------------------------------------------------

# Import Gini index and income data
Gini <- read.csv("./data/Gini_index_worldbank.csv")
income <- read.csv("./data/IncomeGDP.csv")
Gini <- merge(Gini, income, all.x=T)

# Modify the IFR comparison list
country_IFRs <- IFR_comparison %>%
  filter(Region != "Global") %>%
  rename(Country = Region) %>%
  mutate(IFR = IFR/100) %>%
  mutate(Country = str_replace(Country,"New York City", "United States")) %>%
  mutate(Country = str_replace(Country, "UK", "United Kingdom")) %>%
  mutate(bin = NA)

# Categorize IFR estimates into age category bins based on age_classes
agec=c(0,18,45,65,75) #age mins
for (i in 1:nrow(country_IFRs)) {
  country_IFRs$bin[i] = age_classes[max(which(country_IFRs$AgeMidpoint[i]>agec))]
}

# Merge with Gini and income data
country_IFRs_wealth <- aggregate(IFR ~ bin+Country, data=country_IFRs, FUN=mean) %>%
  merge(Gini)

# Create objects for storing income, Gini, and combined
income <- seq(min(country_IFRs_wealth$medianIncome),
              max(country_IFRs_wealth$medianIncome),
              by=50)
Gini <- seq(min(country_IFRs_wealth$Gini),
            max(country_IFRs_wealth$Gini),
            length.out=length(income))
wealth_data <- data.frame(medianIncome=rep(income,length(age_classes)),
                          Gini=rep(Gini,length(age_classes)),
                          bin=rep(age_classes,each=length(income)),
                          Country="")

# Perform regression and find which relationships are significant for plotting
for (i in 1:length(age_classes)) { 
  income_regression <- betareg(IFR~medianIncome,
                               data=country_IFRs_wealth[country_IFRs_wealth$bin==age_classes[i],])
  Gini_regression <- betareg(IFR~Gini,
                             data=country_IFRs_wealth[country_IFRs_wealth$bin==age_classes[i],])
  wealth_data$IFR_inc[wealth_data$bin==age_classes[i]]=
    predict(income_regression, type="response",
            newdata=wealth_data[wealth_data$bin==age_classes[i],])
  wealth_data$IFRse_inc[wealth_data$bin==age_classes[i]]=
    as.data.frame(allEffects(income_regression,
                             xlevels=
                               list(medianIncome=
                                      wealth_data$medianIncome[wealth_data$bin==age_classes[i]])))$medianIncome[,"se"]
  wealth_data$pval_inc[wealth_data$bin==age_classes[i]]=
    summary(income_regression)$coef$mean["medianIncome","Pr(>|z|)"]
  wealth_data$pval_gin[wealth_data$bin==age_classes[i]]=
    summary(Gini_regression)$coef$mean["Gini","Pr(>|z|)"]
  wealth_data$IFR_gin[wealth_data$bin==age_classes[i]]=
    predict(Gini_regression, type="response",
            newdata=wealth_data[wealth_data$bin==age_classes[i],])
  wealth_data$IFRse_gin[wealth_data$bin==age_classes[i]]=
    as.data.frame(allEffects(Gini_regression,
                             xlevels=
                               list(Gini=
                                      wealth_data$Gini[wealth_data$bin==age_classes[i]])))$Gini[,"se"]
}

# Option to see analysis of individual age-category for income, Gini, or both
# i=1: 0-17       i=2: 18-44       i=3: 45-64        i=4: 65-74       i=5: 75+
i=5;summary(betareg(IFR~medianIncome,data=IFRg[IFRg$cat==agecm[i],]));agecm[i]
i=5;summary(betareg(IFR~Gini,data=IFRg[IFRg$cat==agecm[i],]));agecm[i]
i=5;summary(betareg(IFR~medianIncome+Gini,data=IFRg[IFRg$cat==agecm[i],]));agecm[i]

#D rop lowest age bin b/c only 3 data points
country_IFRs_wealth <- filter(country_IFRs_wealth, bin != "0-17" ) 
wealth_data <- filter(wealth_data, bin != "0-17")

# Income plot
ggplot(country_IFRs_wealth, aes(y=IFR, x=medianIncome, label=Country))+
  geom_point()+
  facet_wrap(vars(bin), nrow=2, scales="free_y")+
  theme_few()+
  geom_line(data=wealth_data[wealth_data$pval_inc<0.05,], aes(x=medianIncome,y=IFR_inc))+
  geom_ribbon(data=wealth_data[wealth_data$pval_inc<0.05,],
              aes(x=medianIncome,y=IFR_inc,ymin=IFR_inc-IFRse_inc,ymax=IFR_inc+IFRse_inc),
              alpha=0.5,fill="grey",linetype=0)+
  geom_text_repel(size=4)+
  theme(legend.position = "none")+
  theme(axis.title=element_text(size=25), axis.text=element_text(size=25),
        strip.text.x = element_text(size = 15))+
  scale_y_continuous(labels=scales::percent)+
  scale_x_continuous(breaks=seq(5000,20000,5000),labels=c("5,000","10,000","15,000","20,000"))+
  xlab("Median income ($)")+ 
  ylab("Infection fatality rate")

# Gini index plot
# need to add a b c d to these
# need to decapitalize the labels on the other graphs
ggplot(country_IFRs_wealth,aes(y=IFR, x=Gini, label=Country))+
  geom_point()+
  facet_wrap(vars(bin), nrow=2, scales="free_y")+
  theme_few()+
  geom_line(data=wealth_data[wealth_data$pval_gin<0.05, ], aes(x=Gini, y=IFR_gin))+
  geom_ribbon(data=wealth_data[wealth_data$pval_gin<0.05, ],
              aes(x=Gini,y=IFR_gin,ymin=IFR_gin-IFRse_gin,ymax=IFR_gin+IFRse_gin),
              alpha=0.5,fill="grey",linetype=0)+
  geom_text_repel(size=4)+
  theme(axis.title=element_text(size=20), axis.text=element_text(size=20),
        strip.text.x = element_text(size = 15))+
  scale_y_continuous(labels=scales::percent)+
  xlab("Gini index")+ 
  ylab("Infection fatality rate")
