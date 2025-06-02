## 01_tdr_data_summaries.R: summary statistics of TDR dive data 

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 01 Jun 2025

## --------------------------------------------------------------------------- ##

## load packages 
library(dplyr)
library(here)
library(lubridate)

## read in compiled TDR data ## ---------------------------------------------- ##
tdr <- read.csv(here("data","PcTagTDRDataCompiled.csv"))

# get the total amount of data per individual 
tot_h <- tdr %>%
  group_by(TagID) %>%
  summarise(
    tot_h = sum(DifferenceInTimeNumeric)/60/60
  )
sum(tot_h$tot_h)

# get the records that are less than or equal to 10 m
d10_h <- tdr %>%
  filter(DepthMeter <= 10) %>%
  group_by(TagID) %>%
  summarise(
    d10_h = sum(DifferenceInTimeNumeric)/60/60
  )

# get the records that are less than or equal to 50 m
d50_h <- tdr %>%
  filter(DepthMeter <= 50) %>%
  group_by(TagID) %>%
  summarise(
    d50_h = sum(DifferenceInTimeNumeric)/60/60
  )

# get proportion of time by those two depth categories
prop <- left_join(tot_h, d10_h, by = "TagID") %>%
  left_join(., d50_h) %>%
  mutate(
    prop10 = (d10_h/tot_h)*100,
    prop50 = (d50_h/tot_h)*100
  )

# get mean of proportion in top 10 m
mean(prop$prop10)

# get the number of hours of data by time of day 
tdr %>%
  group_by(TimeOfDay) %>%
  summarise(
    hrs = sum(DurationMinutes)/60
  )

# get dives greater than or equal to 10 m
dives10 <- filter(tdr, DepthMeter >= 10)

# summarize by tag ID
dives10_sum <- dives10 %>%
  group_by(TagID) %>%
  summarise(
    n = n(),
    mean_depth = mean(DepthMeter),
    sd_depth = sd(DepthMeter),
    med_depth = median(DepthMeter),
    max_depth = max(DepthMeter),
    mean_dur = mean(DurationMinutes),
    sd_dur = sd(DurationMinutes),
    med_dur = median(DurationMinutes),
    max_dur = max(DurationMinutes),
    mean_desc = mean(Descent),
    sd_desc = sd(Descent),
    med_desc = median(Descent),
    max_desc = max(Descent),
    mean_asc = mean(abs(Ascent)),
    sd_asc = sd(abs(Ascent)),
    med_asc = median(abs(Ascent)),
    max_asc = max(abs(Ascent))
  )

# get the grand/overall means from individual medians
dives10_all_sum <- dives10_sum %>%
  summarise(
    mean_depth = mean(med_depth),
    sd_depth = sd(med_depth),
    cv_depth = sd(med_depth)/mean_depth,
    mean_dur = mean(med_dur),
    sd_dur = sd(med_dur),
    cv_dur = sd(med_dur)/mean_dur,
    mean_desc = mean(med_desc),
    sd_desc = sd(med_desc),
    cv_desc = sd(med_desc)/mean_desc,
    mean_asc = mean(med_asc),
    sd_asc = sd(med_asc),
    cv_asc = sd(med_asc)/mean_asc
  )

# diel comparisons >50m 
dives50_tod <- tdr %>%
  filter(., DepthMeter >= 50) %>%
  group_by(TagID, TimeOfDay) %>%
  summarise(
    n = n(),
    mean_depth = mean(DepthMeter),
    sd_depth = sd(DepthMeter),
    med_depth = median(DepthMeter),
    max_depth = max(DepthMeter),
    mean_dur = mean(DurationMinutes),
    sd_dur = sd(DurationMinutes),
    med_dur = median(DurationMinutes),
    max_dur = max(DurationMinutes)
  )

# diel comparisons >10m 
dives10_tod <- tdr %>%
  filter(., DepthMeter >= 10) %>%
  filter(., DepthMeter < 230) %>%
  filter(., TagID %in% c("PcTDR02","PcTDR04","PcTDR05")) %>%
  group_by(TagID, TimeOfDay) %>%
  summarise(
    n = n(),
    mean_depth = mean(DepthMeter),
    sd_depth = sd(DepthMeter),
    med_depth = median(DepthMeter),
    max_depth = max(DepthMeter),
    mean_dur = mean(DurationMinutes),
    sd_dur = sd(DurationMinutes),
    med_dur = median(DurationMinutes),
    max_dur = max(DurationMinutes),
    mean_desc = mean(Descent),
    sd_desc = sd(Descent),
    med_desc = median(Descent),
    max_desc = max(Descent),
    mean_asc = mean(abs(Ascent)),
    sd_asc = sd(abs(Ascent)),
    med_asc = median(abs(Ascent)),
    max_asc = max(abs(Ascent))
  ) %>%
  group_by(TimeOfDay) %>%
  summarise(
    mean_depth = mean(med_depth),
    sd_depth = sd(med_depth),
    cv_depth = sd(med_depth)/mean_depth,
    mean_dur = mean(med_dur),
    sd_dur = sd(med_dur),
    cv_dur = sd(med_dur)/mean_dur,
    mean_desc = mean(med_desc),
    sd_desc = sd(med_desc),
    cv_desc = sd(med_desc)/mean_desc,
    mean_asc = mean(med_asc),
    sd_asc = sd(med_asc),
    cv_asc = sd(med_asc)/mean_asc
  )

# get the proportion of time spent in greater than 50 m by time of day for 
# relevant tags 
d50g_h <- tdr %>%
  filter(DepthMeter > 50) %>%
  filter(., DepthMeter < 230) %>%
  filter(., TagID %in% c("PcTDR02","PcTDR04","PcTDR05")) %>%
  group_by(TagID, TimeOfDay) %>%
  summarise(
    d50_h = sum(DifferenceInTimeNumeric)/60/60
  )

# get the total time 
tot_h_tod <- tdr %>%
  filter(., TagID %in% c("PcTDR02","PcTDR04","PcTDR05")) %>%
  group_by(TagID, TimeOfDay) %>%
  summarise(
    tot_h = sum(DifferenceInTimeNumeric)/60/60
  )

# time in >50
prop50_tod <- left_join(tot_h_tod, d50g_h, by = c("TagID","TimeOfDay")) %>%
  mutate(
    prop50 = (d50_h/tot_h)*100
  ) %>%
  filter(., !is.na(prop50)) 

# proportion
prop50_tod_sum <- prop50_tod %>%
  group_by(TimeOfDay) %>%
  summarise(
    mean = mean(prop50),
    sd = sd(prop50),
    cv = sd(prop50)/mean(prop50)
  )


## relative velocity summary stats for relevant tags ## ---------------------- ##

# PcTDR01
pc1 <- read.csv(here("data","PcTDR01_tdr_depth_m_for_r.csv"))
mean(pc1$Velocity, na.rm = T)
sd(pc1$Velocity, na.rm = T)
median(pc1$Velocity, na.rm = T)
max(pc1$Velocity, na.rm = T)

# PcTDR02
pc2 <- read.csv(here("data","PcTDR02_tdr_depth_m_for_r.csv"))
mean(pc2$Velocity, na.rm = T)
sd(pc2$Velocity, na.rm = T)
median(pc2$Velocity, na.rm = T)
max(pc2$Velocity, na.rm = T)

# PcTDR03
pc3 <- read.csv(here("data","PcTDR03_tdr_depth_m_for_r.csv"))
mean(pc3$Velocity, na.rm = T)
sd(pc3$Velocity, na.rm = T)
median(pc3$Velocity, na.rm = T)
max(pc3$Velocity, na.rm = T)

# PcTDR04
# read in both days separately then merge
pc4.1 <- read.csv(here("data","PcTDR04_tdr_depth_m_2004Oct06_for_r.csv"))
summary(pc4.1)
pc4.1$dmy <- "2004-10-06"
pc4.2 <- read.csv(here("data","PcTDR04_tdr_depth_m_2004Oct07_for_r.csv"))
pc4.2$dmy <- "2004-10-07"

# make h:m:s out of numeric time 
pc4.1$hms <- strftime(as.POSIXct(pc4.1$Time * 60 * 60, origin = Sys.Date(), tz = "GMT"), format = "%H:%M:%S") 
pc4.1$datetime <- as.POSIXct(paste0(pc4.1$dmy," ",pc4.1$hms), format = "%Y-%m-%d %H:%M:%S", tz = "Pacific/Honolulu") + hours(7)

pc4.2$hms <- strftime(as.POSIXct(pc4.2$Time * 60 * 60, origin = Sys.Date(), tz = "GMT"), format = "%H:%M:%S") 
pc4.2$datetime <- as.POSIXct(paste0(pc4.2$dmy," ",pc4.2$hms), format = "%Y-%m-%d %H:%M:%S", tz = "Pacific/Honolulu") + hours(7) - days(1)

# combine and add id and lat lon
pc4 <- bind_rows(pc4.1,pc4.2) %>%
  mutate(
    id = "PcTDR04",
    lon = -156.21027,
    lat = 19.60738,
    start_hst = datetime
  )

saveRDS(pc4, here("data","PcTDR04_tdr_depth_m_days_combined_2025Mar18.rds"))

mean(pc4$Velocity, na.rm = T)
sd(pc4$Velocity, na.rm = T)
median(pc4$Velocity, na.rm = T)
max(pc4$Velocity, na.rm = T)

## summary stats for diel period for 02 and 04 that had data over different periods

# add the date to PcTDR02
pc2$dmy <- rep("1999-11-17",nrow(pc2))

# make h:m:s out of numeric time 
pc2$hms <- strftime(as.POSIXct(pc2$Time * 60 * 60, origin = Sys.Date(), tz = "GMT"), format = "%H:%M:%S") 
pc2$datetime <- as.POSIXct(paste0(pc2$dmy," ",pc2$hms), format = "%Y-%m-%d %H:%M:%S", tz = "Pacific/Honolulu")
pc2 <- pc2 %>%
  mutate(datetime = datetime + hours(7))
str(pc2)

# assign a lat/lon (deployment area is fine) for extracting time of day
pc2$lat <- 20.6
pc2$lon <- -156.9
pc2$start_hst <- pc2$datetime
pc2$id <- "PcTDR02"

# combine data for both tags
pc24 <- bind_rows(pc2, pc4)

# extract sun angles and classify time of day 
pc24_sf <- st_as_sf(pc24, coords = c("lon","lat"), crs = 4326) %>%
  mutate(
    sunrise = as.POSIXct(suntools::sunriset(., start_hst, direction = "sunrise", POSIXct.out=T)$time),
    sunset = as.POSIXct(suntools::sunriset(., start_hst, direction = "sunset", POSIXct.out=T)$time),
    solar_noon = as.POSIXct(suntools::solarnoon(., start_hst, POSIXct.out=T)$time),
    civil_dawn = as.POSIXct(suntools::crepuscule(., start_hst, solarDep = 6, direction="dawn", POSIXct.out=T)$time),
    end_dawn = as.POSIXct(suntools::crepuscule(., start_hst, solarDep = -6, direction="dawn", POSIXct.out=T)$time),
    civil_dusk = as.POSIXct(suntools::crepuscule(., start_hst, solarDep = 6, direction="dusk", POSIXct.out=T)$time),
    start_dusk = as.POSIXct(suntools::crepuscule(., start_hst, solarDep = -6, direction="dusk", POSIXct.out=T)$time)
  ) %>%
  mutate(
    tod = case_when(
      start_hst < civil_dawn | start_hst > civil_dusk ~ "night",
      start_hst > end_dawn & start_hst < start_dusk ~ "day",
      start_hst >= civil_dawn & start_hst <= end_dawn ~ "dawn",
      start_hst >= start_dusk & start_hst <= civil_dusk ~ "dusk"
    )
  )
unique(pc24_sf$tod)  

# summarize by id and time of day 
ss_tod <- pc24_sf %>%
  group_by(id, tod) %>%
  summarise(
    n = n(),
    mean_ss = mean(Velocity, na.rm = T),
    median_ss = median(Velocity, na.rm = T),
    max_ss = max(Velocity, na.rm = T),
    cv_ss = sd(Velocity, na.rm = T)/mean_ss
  )

ss_tod_sum <- pc24_sf %>%
  group_by(id, tod) %>%
  summarise(
    n = n(),
    mean_ss = mean(Velocity, na.rm = T),
    median_ss = median(Velocity, na.rm = T)
  ) %>%
  group_by(tod) %>%
  summarise(
    grmean = mean(median_ss),
    cv = sd(median_ss)/grmean
  )
