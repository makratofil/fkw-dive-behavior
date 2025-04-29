## 04_temporal_gamm_models.R: fit GAMM for dive data response variables ~ temporal
## predictor variables 

## in contrast to "spatiotemporal_gamm_models.R", this script uses the whole dive
## dataset (i.e., dive locations not restricted by location error) and only 
## temporal predictors

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 20 Feb 2025

## --------------------------------------------------------------------------- ##

## load packages 
library(dplyr)
library(here)
library(lubridate)
library(ggplot2)
library(mgcv)
library(gratia)
library(gam.hp)

## set seed
set.seed(123)

## read in cleaned and processed behavior log data ## ------------------------ ##
beh <- readRDS(here("pipeline","clean_data_for_analysis",
                    "all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_split_tod_2025Feb18.rds"))

# review to make sure datetimes are already formatted 
str(beh)

# add population column 
beh <- beh %>%
  mutate(
    population = case_when(
      DeployID %in% c("PcTag026","PcTag028","PcTag030","PcTag032","PcTag055",
                      "PcTag074") ~ "MHI",
      DeployID %in% c("PcTag035","PcTag037","PcTag049") ~ "NWHI",
      DeployID %in% c("PcTag090","PcTag092","PcTagP09") ~ "Open-ocean"
    )
  )
beh$population <- as.factor(beh$population)


# filter out dives 
dives <- filter(beh, What == "Dive")
hist(dives$depth_avg50)
hist(dives$dur_mins)
hist(dives$dur_secs)

# select the temporal columns we want to consider for modeling 
dives_temp <- dives %>%
  dplyr::select(
    DeployID, start_hst, end_hst, depth_avg50, dur_secs, dur_mins, population, moon_phase, moon_ill_fraction
  )

# make a continuous version of moon phase
dives_temp$moon_phase_rads <- lunar::lunar.phase(as.Date(dives_temp$start_hst))

# get time in HH:MM:SS
dives_temp$hms <- strftime(dives_temp$start_hst, format= "%H:%M:%S", tz= "Pacific/Honolulu")
dives_temp$hms_n <- as.numeric(hms(dives_temp$hms))

# assess correlation among predictors 
preds <- dives_temp %>%
  dplyr::select(moon_ill_fraction,
                moon_phase_rads,hms_n)
cor(preds) # all under 0.5

# make DeployID a factor leveled by population
dives_temp$DeployID <- factor(dives_temp$DeployID, levels = c("PcTag026","PcTag028","PcTag030","PcTag032","PcTag055",
                                                              "PcTag074","PcTag035","PcTag037","PcTag049",
                                                              "PcTag090","PcTag092","PcTagP09"))

summary(dives_temp$DeployID)

# save RDS file of formatted model data
saveRDS(dives_temp, here("pipeline","gamm_models","temporal_gamms_allpops_model_data_2025Feb20.rds"))


## GAMMs: Depth ## ----------------------------------------------------------- ##

## 1. global smooth: account for non-independence within individuals, and all 
## individuals share a common trend in dive depth ## 
depth_g <- gam(depth_avg50 ~ s(hms_n, bs = "cc", k = 5) +
                 s(moon_phase_rads, bs = "cc", k = 5) +
                 s(DeployID, bs = "re"),
               knots = list(hms_n = c(0,86400), 
                            moon_phase_rads = c(min(dives_temp$moon_phase_rads),
                                                max(dives_temp$moon_phase_rads))),
               family = Gamma(link = "log"),
               data = dives_temp,
               select = T,
               method = "REML")

# review ACF and PACF plots
acf(resid(depth_g, type = "pearson"))
pacf(resid(depth_g, type = "pearson"))

# review outputs
summary(depth_g)
plot(depth_g)
appraise(depth_g)
draw(depth_g, scales = "fixed")

# save model object
saveRDS(depth_g, here("pipeline","gamm_models","temporal_gamms_allpops_depth_g_2025Feb20.rds"))

## 2. global smooth plus tag-specific smooths for time of day (similar penalties) ## 
depth_gs <- gam(depth_avg50 ~ s(hms_n, bs = "cc", k = 5) +
                  s(moon_phase_rads, bs = "cc", k = 5) +
                  s(hms_n, DeployID, bs = "fs", xt=list(bs="cc", k = 5)),
                knots = list(hms_n = c(0,86400), 
                             moon_phase_rads = c(min(dives_temp$moon_phase_rads),
                                                 max(dives_temp$moon_phase_rads))),
                family = Gamma(link = "log"),
                select = T,
                data = dives_temp,
                method = "REML")

# review ACF and PACF plots 
acf(resid(depth_gs, type = "pearson")) 
pacf(resid(depth_gs, type = "pearson"))

# review outputs
summary(depth_gs)
plot(depth_gs)
appraise(depth_gs)
draw(depth_gs, scales = "fixed")

# save model object
saveRDS(depth_gs, here("pipeline","gamm_models","temporal_gamms_allpops_depth_gs_2025Feb20.rds"))

# compare both models
AIC(depth_g, depth_gs) # gs model better fit 

# get the predictor-specific deviance explained
hp_depth_gs <- gam.hp(depth_gs)
hp_depth_gs$hierarchical.partitioning

## GAMMs: Duration ## -------------------------------------------------------- ##

## 1. global smooth: account for non-independence within individuals, and all 
## individuals share a common trend in dive depth ## 
dur_g <- gam(dur_secs ~ s(hms_n, bs = "cc", k = 5) +
               s(moon_phase_rads, bs = "cc", k = 5) +
               s(DeployID, bs = "re"),
             knots = list(hms_n = c(0,86400), 
                          moon_phase_rads = c(min(dives_temp$moon_phase_rads),
                                              max(dives_temp$moon_phase_rads))),
             family = Gamma(link = "log"),
             data = dives_temp,
             select = T,
             method = "REML")

# review ACF and PACF plots
acf(resid(dur_g, type = "pearson"))
pacf(resid(dur_g, type = "pearson"))

# review outputs
summary(dur_g)
plot(dur_g)
appraise(dur_g)
draw(dur_g, scales = "fixed")

# save model object
saveRDS(dur_g, here("pipeline","gamm_models","temporal_gamms_allpops_dur_secs_g_2025Feb20.rds"))

## 2. global smooth plus tag-specific smooths for time of day (similar penalties) ## 
dur_gs <- gam(dur_secs ~ s(hms_n, bs = "cc", k = 5) +
                s(moon_phase_rads, bs = "cc", k = 5) +
                s(hms_n, DeployID, bs = "fs", xt=list(bs="cc", k = 5)),
              knots = list(hms_n = c(0,86400), 
                           moon_phase_rads = c(min(dives_temp$moon_phase_rads),
                                               max(dives_temp$moon_phase_rads))),
              family = Gamma(link = "log"),
              select = T,
              data = dives_temp,
              method = "REML")

# review ACF and PACF plots 
acf(resid(dur_gs, type = "pearson")) 
pacf(resid(dur_gs, type = "pearson"))

# review outputs
summary(dur_gs)
plot(dur_gs)
appraise(dur_gs)
draw(dur_gs, scales = "fixed")

# save model object
saveRDS(dur_gs, here("pipeline","gamm_models","temporal_gamms_allpops_dur_secs_gs_2025Feb20.rds"))

# compare both models
AIC(dur_g, dur_gs) # gs model better fit 

# get the predictor-specific deviance explained
hp_dur_gs <- gam.hp(dur_gs)
hp_dur_gs$hierarchical.partitioning


## GAMMs: Hourly dive rate ## ------------------------------------------------ ##

## read in cleaned and processed behavior log data for HOUR of day ## 
hod <- readRDS(here("pipeline","clean_data_for_analysis",
                    "all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_split_hod_2025Feb18.rds"))

# get column for hour of the day
hod$hod <- hour(hod$start_hst)

# now for the behavior log, get the total number of hours of BL data per hour for 
# each tag 
hod$dmy <- as.Date(hod$start_hst, tz = "Pacific/Honolulu")
hod$dmy_hod <- paste0(hod$dmy, " ", hod$hod)

# sum by hour for all surface and dive records 
hod_sum <- hod %>%
  group_by(DeployID, dmy_hod) %>%
  summarise(
    hrs = sum(dur_hrs)
  )

# filter out dives 
dives <- filter(hod, What == "Dive")

# get the number of dives per hour per day for each tag 
dives_hod <- dives %>%
  group_by(DeployID, dmy_hod) %>%
  summarise(
    n_dives = n(),
    avg_depth = mean(depth_avg50),
    max_depth = max(depth_avg50),
    avg_dur = mean(dur_mins)
  )

# combine the two; any records without a dive can be assigned 0 for no dives
rate <- left_join(hod_sum, dives_hod, by = c("DeployID","dmy_hod")) %>%
  mutate(
    n_dives = ifelse(is.na(n_dives),0, n_dives),
    avg_depth = ifelse(is.na(avg_depth),0, avg_depth),
    max_depth = ifelse(is.na(max_depth),0, max_depth),
    avg_dur = ifelse(is.na(avg_dur),0,avg_dur),
    dph = n_dives/hrs
  )

# remove hours that have <75% of data 
rate_sub <- filter(rate, hrs >= 0.75)
hist(rate_sub$dph)
hist(rate_sub$avg_depth)

# get the covariates of interest for the dmy 
hod_covs <- hod %>%
  group_by(DeployID, dmy_hod) %>%
  slice(1) %>%
  select(DeployID, dmy_hod, dmy, hod)

# add them to the rates
rate_temp <- left_join(rate_sub, hod_covs, by = c("DeployID","dmy_hod"))

# quick summary by ID
rate_temp_sum <- rate_temp %>%
  group_by(DeployID, hod) %>%
  tally()

# make a continuous version of moon phase
rate_temp$moon_phase_rads <- lunar::lunar.phase(rate_temp$dmy)

# assess correlation among predictors 
preds <- rate_temp %>%
  ungroup() %>%
  dplyr::select(hod,
                moon_phase_rads)
cor(preds) # all under 0.5

# make DeployID a factor leveled by population
rate_temp$DeployID <- factor(rate_temp$DeployID, levels = c("PcTag026","PcTag028","PcTag030","PcTag032","PcTag055",
                                                            "PcTag074","PcTag035","PcTag037","PcTag049",
                                                            "PcTag090","PcTag092","PcTagP09"))

summary(rate_temp$DeployID)

# save data 
saveRDS(rate_temp, here("pipeline","gamm_models","temporal_gamms_rate_model_data_2025Feb20.rds"))


## GAMM: hourly rate ~ temporal variables ## -------------------------------- ##

## 1. global smooths and random intercept for ID. fit with negative binomial ## 
rate_g <- gam(n_dives ~ s(hod, bs = "cc", k = 5) +
                s(moon_phase_rads, bs = "cc", k = 5) +
                s(DeployID, bs = "re") +
                offset(log(hrs)),
              knots = list(hod = c(0,23), 
                           moon_phase_rads = c(min(rate_temp$moon_phase_rads),
                                               max(rate_temp$moon_phase_rads))),
              family = nb(link = "log"),
              select = T,
              data = rate_temp,
              method = "REML"
)
# review ACF and PACF plots 
acf(resid(rate_g, type = "pearson")) 
pacf(resid(rate_g, type = "pearson"))

# review
summary(rate_g)
plot(rate_g)
draw(rate_g, scales = "fixed")
appraise(rate_g)
draw(rootogram(rate_g)) # recommended to look at this instead of standing histogram Kleiber & Zeilis 2016

# save model object 
saveRDS(rate_g, here("pipeline","gamm_models","temporal_gamms_allpops_rate_g_2024Dec23.rds"))

## 2. factor-level smooths for hour of the day ## 
rate_gs <- gam(n_dives ~ s(hod, bs = "cc", k = 5) +
                 s(moon_phase_rads, bs = "cc", k = 5) +
                 s(hod, DeployID, bs = "fs", xt=list(bs="cc", k = 5)) +
                 offset(log(hrs)),
               knots = list(hod = c(0,23), 
                            moon_phase_rads = c(min(rate_temp$moon_phase_rads),
                                                max(rate_temp$moon_phase_rads))),
               family = nb(link = "log"),
               select = T,
               data = rate_temp,
               method = "REML"
)

# review ACF and PACF plots 
acf(resid(rate_gs, type = "pearson")) 
pacf(resid(rate_gs, type = "pearson"))

# review
summary(rate_gs)
draw(rate_gs, scales = "fixed")
appraise(rate_gs)
draw(rootogram(rate_gs)) # recommended to look at this instead of standing histogram Kleiber & Zeilis 2016

# save model object 
saveRDS(rate_gs, here("pipeline","gamm_models","temporal_gamms_allpops_rate_gs_2025Feb20.rds"))

# compare models 
AIC(rate_g, rate_gs)

# note: gam.hp doesn't have functionality to deal with models that have offsets