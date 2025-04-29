## 03_spatiotemporal_gamm_models.R: fit GAMM for dive data response variables ~ 
## spatial and temporal predictor variables

## all populations, models only (prediction plots in separate script)

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 21 Feb 2025

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
summary(beh)

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

# filter out dives that have estimated error > 4km
dives <- beh %>%
  filter(What == "Dive") %>%
  mutate(
    avg_se = (se.mu.x + se.mu.y)/2,
    dive_num = seq(1,n(), by = 1)
  ) %>%
  filter(., avg_se < 4000) # removes about 300 dives

# save this (for seafloor dives analysis)
saveRDS(dives, here("pipeline","fkw_dive_ids_4kmerror_for_seafloor_analysis.rds"))

# quick review
summary(dives)
hist(dives$depth_avg50)
hist(dives$dur_mins)
hist(dives$dur_secs)

# look at the number of dives per ID after filtering 
dives %>%
  group_by(DeployID) %>%
  tally()

# rename slope variable (use the lower resolution/GEBCO since it covers all 
# individuals' range)
dives <- dives %>%
  mutate(
    sf_slope = slopeL
  )

# remove dives with chla values that have been flagged with "1" - this flag 
# indicates cells overlapping with land. remove mld values with NA (similar reasons)
unique(dives$chla_flag)
dives_sub <- filter(dives, chla_flag != 1) %>%
  filter(., !is.na(mld)) # removes 14 dives

# take the log10 of chlorophyll-a
dives_sub$logchla <- log10(dives_sub$chla)

# select the columns we want to consider for modeling 
dives_spat <- dives_sub %>%
  dplyr::select(
    DeployID, start_hst, end_hst, depth_avg50, dur_mins, dur_secs, population, moon_phase,
    tod, mld, logchla, sf_slope)

# make a continuous version of moon phase
dives_spat$moon_phase_rads <- lunar::lunar.phase(as.Date(dives_spat$start_hst))

# get time in HH:MM:SS
dives_spat$hms <- strftime(dives_spat$start_hst, format= "%H:%M:%S", tz= "Pacific/Honolulu")
dives_spat$hms_n <- as.numeric(hms(dives_spat$hms))

# assess correlation among predictors 
preds <- dives_spat %>%
  dplyr::select(sf_slope,
                moon_phase_rads,hms_n, logchla,mld) 
cor(preds) # all under 0.5

## simple univariate plots ## ------------------------------------------------ ##

# dive depth vs slope 
ggplot(dives_spat, aes(x = sf_slope, y = depth_avg50)) +
  geom_point() 

# dive depth vs chla 
ggplot(dives_spat, aes(x = logchla, y = depth_avg50)) +
  geom_point()

# dive depth vs mld
ggplot(dives_spat, aes(x = mld, y = depth_avg50)) +
  geom_point() 

# dive depth vs time of day 
ggplot(dives_spat, aes(x = hms_n, y = depth_avg50)) +
  geom_point() +
  facet_wrap(~DeployID, scales = "free")

# dive depth vs moon phase  
ggplot(dives_spat, aes(x = moon_phase_rads, y = depth_avg50)) +
  geom_point() 

## model data prep ## -------------------------------------------------------- ##

# rename 
dives_pop <- dives_spat 

# make ID a factor for random effect
dives_pop$DeployID <- factor(dives_pop$DeployID, levels = c("PcTag026","PcTag028","PcTag030","PcTag032","PcTag055",
                                                            "PcTag074","PcTag035","PcTag037","PcTag049","PcTag090",
                                                            "PcTag092","PcTagP09"))

summary(dives_pop$DeployID) # check

# save RDS file of formatted model data
saveRDS(dives_pop, here("pipeline","gamm_models","spatiotemporal_gamms_allpops_model_ff_data_2025Feb21.rds"))

## dive depth ## ------------------------------------------------------------- ##

## 1. global smooth: account for non-independence within individuals, and all 
## individuals share a common trend in dive depth ## 
depth_pop_g <- gam(depth_avg50 ~ s(hms_n, bs = "cc", k = 5) +
                     s(moon_phase_rads, bs = "cc", k = 5) +
                     s(sf_slope, bs = "tp", k = 5) +
                     s(logchla, bs = "tp", k = 5) +
                     s(mld, bs = "tp", k = 5) +
                     s(DeployID, bs = "re"),
                   knots = list(hms_n = c(0,86400), 
                                moon_phase_rads = c(min(dives_pop$moon_phase_rads),
                                                    max(dives_pop$moon_phase_rads))),
                   family = Gamma(link = "log"),
                   data = dives_pop,
                   select = T,
                   method = "REML")

# review ACF and PACF plots
acf(resid(depth_pop_g, type = "pearson"))
pacf(resid(depth_pop_g, type = "pearson"))

# review outputs
summary(depth_pop_g)
plot(depth_pop_g)
appraise(depth_pop_g)
draw(depth_pop_g, scales = "fixed")

# save model object
saveRDS(depth_pop_g, here("pipeline","gamm_models","spatiotemporal_gamms_allpops_ff_depth_pop_g_2025Feb21.rds"))

## 2. global smooth plus tag-specific smooths for time of day (similar penalties) ## 
depth_pop_gs <- gam(depth_avg50 ~ s(hms_n, bs = "cc", k = 5) +
                      s(moon_phase_rads, bs = "cc", k = 5) +
                      s(sf_slope, bs = "tp", k = 5) +
                      s(logchla, bs = "tp", k = 5) +
                      s(mld, bs = "tp", k = 5) +
                      s(hms_n, DeployID, bs = "fs", xt=list(bs="cc", k = 5)),
                    knots = list(hms_n = c(0,86400), 
                                 moon_phase_rads = c(min(dives_pop$moon_phase_rads),
                                                     max(dives_pop$moon_phase_rads))),
                    family = Gamma(link = "log"),
                    data = dives_pop,
                    select = T,
                    method = "REML")

# review ACF and PACF plots
acf(resid(depth_pop_gs, type = "pearson"))
pacf(resid(depth_pop_gs, type = "pearson"))

# review outputs
summary(depth_pop_gs)
plot(depth_pop_gs)
appraise(depth_pop_gs)
draw(depth_pop_gs, scales = "fixed")

# save model object
saveRDS(depth_pop_gs, here("pipeline","gamm_models","spatiotemporal_gamms_allpops_ff_depth_pop_gs_2025Feb21.rds"))

# compare both models
AIC(depth_pop_g, depth_pop_gs) # gs model better fit 

# get the predictor-specific deviance explained
hp_depth_pop_gs <- gam.hp(depth_pop_gs)
hp_depth_pop_gs$hierarchical.partitioning

## dive duration ## ---------------------------------------------------------- ##
dur_pop_g <- gam(dur_secs ~ s(hms_n, bs = "cc", k = 5) +
                   s(moon_phase_rads, bs = "cc", k = 5) +
                   s(sf_slope, bs = "tp", k = 5) +
                   s(logchla, bs = "tp", k = 5) +
                   s(mld, bs = "tp", k = 5) +
                   s(DeployID, bs = "re"),
                 knots = list(hms_n = c(0,86400), 
                              moon_phase_rads = c(min(dives_pop$moon_phase_rads),
                                                  max(dives_pop$moon_phase_rads))),
                 family = Gamma(link = "log"),
                 data = dives_pop,
                 select = T,
                 method = "REML")

# review ACF and PACF plots
acf(resid(dur_pop_g, type = "pearson"))
pacf(resid(dur_pop_g, type = "pearson"))

# review outputs
summary(dur_pop_g)
plot(dur_pop_g)
appraise(dur_pop_g)
draw(dur_pop_g, scales = "fixed")

# save model object
saveRDS(dur_pop_g, here("pipeline","gamm_models","spatiotemporal_gamms_allpops_ff_duration_pop_g_2025Feb21.rds"))


## 2. global smooth plus tag-specific smooths for time of day (similar penalties) ## 
dur_pop_gs <- gam(dur_secs ~ s(hms_n, bs = "cc", k = 5) +
                    s(moon_phase_rads, bs = "cc", k = 5) +
                    s(sf_slope, bs = "tp", k = 5) +
                    s(logchla, bs = "tp", k = 5) +
                    s(mld, bs = "tp", k = 5) +
                    s(hms_n, DeployID, bs = "fs", xt=list(bs="cc", k = 5)),
                  knots = list(hms_n = c(0,86400), 
                               moon_phase_rads = c(min(dives_pop$moon_phase_rads),
                                                   max(dives_pop$moon_phase_rads))),
                  family = Gamma(link = "log"),
                  data = dives_pop,
                  select = T,
                  method = "REML")


# review outputs
summary(dur_pop_gs)
plot(dur_pop_gs)
appraise(dur_pop_gs)
draw(dur_pop_gs, scales = "fixed")

# save model object
saveRDS(dur_pop_gs, here("pipeline","gamm_models","spatiotemporal_gamms_allpops_ff_duration_pop_gs_2025Feb21.rds"))

# compare both models
AIC(dur_pop_g, dur_pop_gs) # gs model better fit 

# get the predictor-specific deviance explained
hp_dur_pop_gs <- gam.hp(dur_pop_gs)
hp_dur_pop_gs$hierarchical.partitioning