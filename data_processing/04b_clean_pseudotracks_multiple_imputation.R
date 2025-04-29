## clean_pseudotracks_multiple_imputation.R: clean imputed pseudotracks dataset, that has
## already been geoprocessed, for further analysis on dive depth relative
## to seafloor depth

## 04b: multiple imputation

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 23 Feb 2025

## --------------------------------------------------------------------------- ##

## load packages 
library(dplyr)
library(lubridate)

## read in compiled pseudotrack file and format for splitting surfaces ## ---- ##
beh <- readRDS(here("pipeline","geoprocessed",
                    "all_behavlog_pseudotracks_20imputes_rerouted20mIso_geoprocessed_seafloor_2025Feb23.rds"))

# first: for analysis, we will consider dives as 50 meters or greater and 2 minutes
# or longer to maintain consistency across programming regimes. therefore, anything
# less than these thresholds will be considered a surface period 

# update the duration values based on timestamps 
beh <- beh %>%
  mutate(
    end_utc = as.POSIXct(end_utc, tz = "UTC"),
    dur_secs = as.integer(as.duration(start_utc %--% end_utc)),
    dur_mins = dur_secs/60,
    dur_hrs = dur_mins/60,
    dur_days = dur_hrs/24
  )

# get the average depth from the minimum and maximum estimates
beh$depth_avg <- (beh$DepthMin + beh$DepthMax) / 2

# replace depths that are less than 50 with NA and assign as surface periods. 
# make a separate depth column first to check that it worked correctly
beh <- beh %>%
  mutate(depth_avg50 = replace (depth_avg, depth_avg<50, NA)) %>%
  mutate(What = ifelse (is.na(depth_avg50), "Surface", "Dive")) 

# now we need to assign depths of 50m+ but less than 2mins as surface periods.
# first filter out those that meet this threshold, so that you can check that 
# the code for applying the threshold works correctly. for the tags programmed 
# with a 120 second (2min) threshold (PcTag026-037), get the minimum duration 
# in mins to use as the threshold, as this may vary by a very small amount.
min_dur_sum <- beh %>%
  filter(What == "Dive") %>%
  group_by(DeployID) %>%
  summarise(
    min_dur = min(dur_mins)
  ) # use a threshold of 1.9 mins

chk <- filter(beh, depth_avg50 > 50 & dur_mins < 1.9)
summary(chk$depth_avg50, na.rm = T)
summary(chk$dur_mins, na.rm = T)

# now apply the code for the last threshold 
beh2 <- beh %>%
  mutate(depth_avg50 = ifelse (depth_avg50 > 50 & dur_mins < 1.9, NA, depth_avg50)) %>%
  mutate(What = 
           ifelse (is.na(depth_avg50), "Surface", "Dive")) 

# filter out one of the tags to see if it worked 
pc49 <- filter(beh2, DeployID == "PcTag049")

# seems to have worked, so clean up the data a little bit: first remove the 
# original depth and duration columns as to not mix them up further along in the
# analysis 

# for imputed pseudotracks 
beh_clean <- beh2 %>%
  dplyr::select(DeployID, Ptt, crawl_timestamp, trk_id, Start, End, What, Shape, depth_avg50,
                dur_secs, dur_mins, dur_hrs, dur_days, lon, lat, 
                start_utc, end_utc,
                start_hst, depthH, slopeH, aspectH, 
                depthM, slopeM, aspectM, depthL, slopeL, aspectL, sunrise, sunset,
                solar_noon, civil_dawn, end_dawn, civil_dusk, start_dusk, sun_azimuth,
                sun_altitude, moon_azimuth, moon_altitude, moon_ill_fraction, moon_phase,
                tod)

## for analysis of dive depth relative to seafloor depth: we can remove surface
## periods, because we are only looking at dives ## -------------------------- ##
dives <- filter(beh_clean, What == "Dive")

## lastly: remove dives that don't have high resolution (50m MHI or 60m Falkor
## for NWHI) ----------------------------------------------------------------- ##
hires <- filter(dives, !is.na(depthH) | !is.na(depthM))
unique(hires$DeployID)

## save the file ## ---------------------------------------------------------- ##
write.csv(hires, here("pipeline","clean_data_for_analysis",
                      "sub_behavlog_pseudotracks_20imputes_for_seafloor_analysis_2025Feb23.csv"), row.names = F)
saveRDS(hires, here("pipeline","clean_data_for_analysis",
                    "sub_behavlog_pseudotracks_20imputes_for_seafloor_analysis_2025Feb23.rds"))
