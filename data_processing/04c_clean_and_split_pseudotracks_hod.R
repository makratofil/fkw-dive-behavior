## clean_and_split_pseudotracks_hod.R: split surface periods by hour of the day
## in pseudotrack files

## 04c: split by hour of the day, instead of time of day (diel period)

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 10 Oct 2025

## --------------------------------------------------------------------------- ##

## load packages 
library(dplyr)
library(here)
library(lubridate)
library(purrr)

#install.packages("remotes")
#remotes::install_github("NicChr/timeplyr")
library(timeplyr)

## read in compiled pseudotrack file and format for splitting surfaces ## ---- ##
beh <- readRDS(here("pipeline","geoprocessed","all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_2025Oct17.rds"))

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
                sunrise, sunset, solar_noon, civil_dawn, end_dawn, civil_dusk, start_dusk,
                depthH, slopeH, depthM, slopeM, depthL, slopeL, moon_ill_fraction, 
                moon_phase, tod, chlad, chlad_unc, chlad_flag,chla30d, chla30d_unc,
                chla30d_flag, pp, pp_unc, pp_flag, sst, sst_sd, ssh, ssh_sd, 
                u_cur, v_cur, mld)

# get the number of dive and surface periods not split by time of day
beh_clean %>%
  group_by(What) %>%
  tally()

# format to combine with splitted surface data 
beh_clean$end_hst <- as.POSIXct(format(beh_clean$end_utc, tz = "Pacific/Honolulu"), tz = "Pacific/Honolulu")

## now do the splits ## ----------------------------------------------------- ##

# since no dive is going to last that long, splitting surface periods only is fine
dives <- filter(beh_clean, What == "Dive")
surfs <- filter(beh_clean, What == "Surface")

# get the hour day of the start and end of the surface period, and calculate
# the number of hours different between start and end
surfs <- surfs %>%
  mutate(
    hr_day_st = hour(start_hst),
    hr_day_end = hour(end_hst),
    hr_day_diff = hr_day_end - hr_day_st
  )


# create empty dataframe with all of the surface columns
surfs <- as.data.frame(surfs)
bh <- surfs[0,]

# loop through surface periods and determine whether the surface period includes
# more than 1 hour day bin. If so, split up the surface period by the hour-day bins. 
for(i in 1:nrow(surfs)) {
  cat("\r", i, "                     ")
  
  # for testing
  #i = 3
  
  # bind empty dataframe with ith row of surface data
  bh = bind_rows(bh, surfs[i, ])
  
  start <- surfs[i, "start_hst"]
  end <- surfs[i, "end_hst"]
  hr_day_diff <- surfs[i, "hr_day_diff"]
  id <- surfs[i, "DeployID"]
  tod <- surfs[i, "tod"]
  segid4 <- surfs[i, "seg_id4"]
  segid6 <- surfs[i, "seg_id6"]
  segid8 <- surfs[i, "seg_id8"]
  
  # set deploy time zone
  if(id %in% c("PcTag035","PcTag037")) {
    tag_tz = "Pacific/Honolulu"
  } else{
    tag_tz = "UTC"
  }
  
  # determine whether the surface period crosses more than 1 hour (clock hour)
  if(abs(hr_day_diff) > 0) {
    
    # get the hours (rounded) that will be used to split intervals
    times_split <- time_seq(ceiling_date(start, unit = "hour"), end, time_by = "hour")
    
    # create object with start and end times to split surface intervals with
    split_dates <- c(start, times_split, end)
    
    # get intervals of splits 
    splits <- as.data.frame(int_diff(split_dates))
    splits$start <- int_start(splits$`int_diff(split_dates)`)
    splits$end <- int_end(splits$`int_diff(split_dates)`)
    
    # shape up 
    splits_df <- splits %>%
      rename(
        start_hst = start,
        end_hst = end
      ) %>%
      mutate(
        #Start = as.POSIXct(format(start_hst, tz = tag_tz), tz = tag_tz),
        #End = as.POSIXct(format(end_hst, tz = tag_tz), tz = tag_tz),
        DeployID = id,
        tod = tod,
        seg_id4 = segid4,
        seg_id6 = segid6,
        seg_id8 = segid8,
        What = "Surface"
      )
    
    #bh[nrow(bh), "End"] <- as.POSIXct(format(splits_df[1, "end_hst"], tz = tag_tz), tz = tag_tz)
    bh[nrow(bh), "start_hst"] <- splits_df[1,"start_hst"]
    bh[nrow(bh), "end_hst"] <- splits_df[1,"end_hst"]
    
    # now append the rest of the splits
    bh <- bind_rows(bh, splits_df[2:nrow(splits_df),]) %>%
      mutate(
        dur_secs = as.integer(as.duration(start_hst %--% end_hst)),
        dur_mins = dur_secs/60,
        dur_hrs = dur_mins/60,
        dur_days = dur_hrs/24,
        hr_day = hour(start_hst)
      ) %>%
      dplyr::select(-`int_diff(split_dates)`)
    
  } 
  
  
}

# assign to split surface object 
surfs_split <- bh

# format the utc and hst columns of the splitted surface periods 
surfs_split$start_utc <- as.POSIXct(format(surfs_split$start_hst, tz = "UTC"), tz = "UTC")
surfs_split$end_utc <- as.POSIXct(format(surfs_split$end_hst, tz = "UTC"), tz = "UTC")

# remove start and end times from dives (for merging the two datasets)
dives2 <- dives %>% 
  dplyr::select(-c(Start,End))

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
  dplyr::select(-c( Start, End)) %>%
  dplyr::select(DeployID, Ptt, start_utc, end_utc, What, depth_avg50, dur_secs, dur_mins, dur_hrs, dur_days,
                lon, lat, se.mu.x, se.mu.y, seg_id4, seg_id6, seg_id8, start_hst, end_hst, sunrise, sunset,
                civil_dawn, end_dawn, civil_dusk, start_dusk, 
                moon_ill_fraction, moon_phase, tod, 
                depthH, slopeH,depthM, slopeM, depthL, slopeL, Shape,chlad, chlad_unc, chlad_flag,chla30d, chla30d_unc,
                chla30d_flag, pp, pp_unc, pp_flag, sst, sst_sd, ssh, ssh_sd, 
                u_cur, v_cur, mld)

# review
summary(all)

## save the file ## ---------------------------------------------------------- ##
write.csv(all, here("pipeline","clean_data_for_analysis",
                    "all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_split_hod_2025Oct17.csv"), row.names = F)
saveRDS(all, here("pipeline","clean_data_for_analysis",
                  "all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_split_hod_2025Oct17.rds"))
