## clean_and_split_pseudotrack.R: standardize dive/surface periods based on programming
## settings, and split surface periods by time of day (dawn, day, dusk or night) 
## in pseudotrack files

## 04a: regular/single pseudotracks

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 18 Feb 2025

## --------------------------------------------------------------------------- ##

## load packages 
library(dplyr)
library(here)
library(lubridate)
library(purrr)

## read in compiled pseudotrack file and format for splitting surfaces ## ---- ##
beh <- readRDS(here("pipeline","geoprocessed",
                    "all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_2025Feb18.rds"))

# first: for analysis, we will consider dives as 50 meters or greater and 2 minutes
# or longer to maintain consistency across programming regimes. therefore, anything
# less than these thresholds will be considered a surface period 

# get the mean value between the two estimates for dive depth and duration. 
beh <- beh %>%
  mutate(
    depth_avg = if_else(What == "Dive", (DepthMin + DepthMax) / 2, NA_real_),
    dur_secs = (DurationMin + DurationMax) / 2,
    dur_mins = dur_secs/60,
    dur_hrs = dur_mins/60,
    dur_days = dur_hrs/24,
    end_utc = as.POSIXct(end_utc, tz = "UTC")
  )

# replace depths that are less than 50 with NA and assign as surface periods. 
# make a separate depth column first to check that it worked correctly
beh <- beh %>%
  mutate(depth_avg50 = replace (depth_avg, depth_avg < 50, NA)) %>%
  mutate(What = ifelse(is.na(depth_avg50), "Surface", "Dive")) 

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

# quick check; only tags that were programmed to have shorter duration thresholds
# should be included in the 'chk' object
chk <- filter(beh, depth_avg50 > 50 & dur_mins < 1.9)
summary(chk$depth_avg50, na.rm = T)
summary(chk$dur_mins, na.rm = T)

# now apply the code for the last threshold 
beh2 <- beh %>%
  mutate(depth_avg50 = ifelse (depth_avg50 > 50 & dur_mins < 1.9, NA, depth_avg50)) %>%
  mutate(What = 
           ifelse (is.na(depth_avg50), "Surface", "Dive")) 

# filter out one of the tags in 'chk' above to see if it worked 
pc49 <- filter(beh2, DeployID == "PcTag049")

# seems to have worked, so clean up the data a little bit: first remove the 
# original depth and duration columns as to not mix them up further along in the
# analysis 
beh_clean <- beh2 %>%
  dplyr::select(DeployID, Ptt, Start, End, What, Shape, depth_avg50,
                dur_secs, dur_mins, dur_hrs, dur_days, lon, lat, se.mu.x,
                se.mu.y, seg_id4, seg_id6, seg_id8, start_utc, end_utc, start_hst,
                depthH, slopeH, aspectH, depthM, slopeM, aspectM, depthL, slopeL, aspectL,
                sunrise, sunset, solar_noon, civil_dawn, end_dawn, civil_dusk, start_dusk, 
                sun_azimuth, sun_altitude, moon_azimuth, moon_altitude, moon_ill_fraction, 
                moon_phase, tod, chla, chla_unc, chla_flag, mld)

# get the number of dive and surface periods not split by time of day
beh_clean %>%
  group_by(What) %>%
  tally()

## now do the splits ## ----------------------------------------------------- ##

# source surface splitting helper functions
source(here("code","data_processing","surface_splitting_helper_functions.R"))

# need to add this for the function to work 
beh_clean$datetime_hst <- beh_clean$start_hst

# get surface periods 
surfaces <- dplyr::filter(beh_clean, What == "Surface")  %>%
  mutate(
    durationDays = as.duration(start_utc %--% end_utc)/ddays(1)
  )

# see how long the longest surface period is, so that you can verify when you review 
# the splitted surface periods 
summary(surfaces$dur_days)
summary(surfaces$durationDays)

# make a dataframe of time zones that each tag was programmed in (will need this
# for surface splitting functions)
tz_df <- data.frame(
  DeployID = c("PcTag026","PcTag028","PcTag030","PcTag032","PcTag035","PcTag037",
               "PcTag049","PcTag055","PcTag074","PcTag090","PcTag092","PcTagP09"),
  TZ = c(rep("UTC",4), rep("Pacific/Honolulu",2), rep("UTC",6))
)

# add timezone info to the surfs data 
surfaces <- left_join(surfaces, tz_df, by = "DeployID")

# create duplicate DeployID to nest data by 
surfaces$depid <- surfaces$DeployID

# nest dataframe 
surfs_nest <- surfaces %>%
  group_by(depid) %>%
  tidyr::nest()

# create helper functions to create columns in nested dataframe that will then
# be used to apply split surfaces functions to each behavior log 
surfs_nest$TZ <- c(rep("UTC",4), rep("Pacific/Honolulu",2), rep("UTC",6))

format_time <- function(data, TZ){
  
  data$Start <- as.POSIXct(data$Start, tz = TZ)
  data$End <- as.POSIXct(data$End, tz = TZ)
  return(data)
  
}

create_empty <- function(data){
  
  empty <- data[0,]
  return(empty)
}

# apply the functions and split surfaces 
surfs_nest <- surfs_nest %>%
  mutate(
    data2 = purrr::map2(data, TZ, ~format_time(.x, .y)),
    bh = purrr::map(data2, ~create_empty(.x)),
    splits = purrr::map2(data2, bh, ~split_surfs_tod(surfs = .x, bh = .y))
  )

# unnest the now split surface periods 
surfs_split <- surfs_nest %>%
  dplyr::select(splits) %>%
  tidyr::unnest()

# review structure
str(surfs_split)

# check time zone of HST tags
pc35_orig <- filter(surfaces, DeployID == "PcTag035")
pc35_chk <- filter(surfs_split, DeployID == "PcTag035")
tz(pc35_chk$Start)

pc37_orig <- filter(surfaces, DeployID == "PcTag037")
pc37_chk <- filter(surfs_split, DeployID == "PcTag037")
tz(pc35_chk$Start)

# looks like the Start and End times were formatted to UTC. They are correct, 
# and the splitting was done correctly, but will only keep columns in the final
# analytical file with UTC and HST timestamps as to avoid confusion moving 
# forward (i.e., remove Start/End columns)

## merge split surface periods with dives ## --------------------------------- ##

# first get the dives
dives <- filter(beh_clean, What == "Dive")
str(dives)

# format to combine with splitted surface data 
dives$start_hst <- dives$datetime_hst
dives$end_hst <- as.POSIXct(format(dives$end_utc, tz = "Pacific/Honolulu"), tz = "Pacific/Honolulu")

# remove start and end times from dives
dives2 <- dives %>% 
  dplyr::select(-c(Start,End))

# format the utc and hst columns of the splitted surface periods 
surfs_split$start_utc <- as.POSIXct(format(surfs_split$start_hst, tz = "UTC"), tz = "UTC")
surfs_split$end_utc <- as.POSIXct(format(surfs_split$end_hst, tz = "UTC"), tz = "UTC")

# merge with splitted data, fill in column for duration in mins, hours, days,
# and remove columns we don't need or want 
colnames(dives2)
colnames(surfs_split)

all <- bind_rows(dives2, surfs_split) %>%
  arrange(DeployID, start_hst) %>%
  mutate(
    dur_secs = as.integer(as.duration(start_hst %--% end_hst)),
    dur_mins = dur_secs/60,
    dur_hrs = dur_mins/60,
    dur_days = dur_hrs/24
  ) %>%
  dplyr::select(-c(durationDays, `int_diff(dates)`, Start, End, datetime_hst, depid, TZ)) %>%
  dplyr::select(DeployID, Ptt, start_utc, end_utc, What, depth_avg50, dur_secs, dur_mins, dur_hrs, dur_days,
         lon, lat, se.mu.x, se.mu.y, seg_id4, seg_id6, seg_id8, start_hst, end_hst, sunrise, sunset,
         civil_dawn, end_dawn, civil_dusk, start_dusk, sun_azimuth, sun_altitude, moon_azimuth, 
         moon_altitude, moon_ill_fraction, moon_phase, tod, 
         depthH, slopeH,aspectH, depthM, slopeM, aspectM, depthL, slopeL, aspectL, Shape,
         chla, chla_unc, chla_flag, mld)

# review
summary(all)

## save the file ## ---------------------------------------------------------- ##
write.csv(all, here("pipeline","clean_data_for_analysis",
                        "all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_split_tod_2025Feb18.csv"), row.names = F)
saveRDS(all, here("pipeline","clean_data_for_analysis",
                      "all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_split_tod_2025Feb18.rds"))
