READ ME

This is the code associated with Rickards & Kilpatrick 2022.

If you wish to replcate the results, run the files in the following order:
1. draw_ifrs.R
- This creates the seroprevalence draws and runs the statistical test that generates the IFRs
2. postprocess.R
- This processes the IFR draws into result tables, and provides graphs of the delay distributions, epidata, and timeline.


Additional scripts:

"edit_serodata.R" reformats the serodata from Rosenberg et al. 2020 into the age groups used by the New York City Department of Health and Mental Hygiene

"extend_dates.R" prepends age-stratified case and death data based on non-stratified case and death counts, for the first two weeks of the pandemic

"ifr_betabinomial.stan" contains the STAN file needed to carry out the Bayesian inference

"res_file.rds" extracts the IFR from the "stanfit" file

"stanfit_file.rds" contains the results from sampling from the Bayesian model

"utils.R" contains additional functions needed in "draw_ifrs.R" and "postprocess.R"


Data objects and miscellaneous:

"age_stratified_IFRs_[suffix].rds" -  contains all of the IFRs estimations, not yet processed, after it's gone through the stan model

"data_setup.rda" - carries data over from "draw_ifrs" to "postprocess"


Files in "data" folder:

"Gini_index_worldbank.csv" contains the Gini index for the countries of interest

"IFR_comparison.csv" contains the IFRs for the countries of interest

"IncomeGDP.csv" contains the income and GDP information for the countries of interest

"NYC_pop_age_total.csv" contains the population for each age class of interest

"NYCHealth_cumul_data_extd.csv" contains the edited and extended case and death data, cumulative for each age class up to the day

"NYCHealth_cumul_data.csv" contains the unedited and unextended case and death data, cumulative for each age class up to the day

"NYCHealth_daily_data.csv" contains daily case and death data

"serodata_edited.csv" contains serodata in the age classes reported by the NYC Department of Health and Mental Hygiene

"serodata.csv" contains serodata in the age classes originally reported by Rosenberg et al. 2020

"serop_age_spain.csv" contains seroprevalence data from Pollan et al. 2020, which is used here to calibrate the serodata for the 0-17 age class


Files in "results" folder contain various tables and graphs created by "postprocess.R"


