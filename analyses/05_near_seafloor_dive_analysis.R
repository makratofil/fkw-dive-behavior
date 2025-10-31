## 05_near_seafloor_dive_analysis.R: assess dive depth ~ seafloor depth for 
## insular animals 

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 19 Oct 2025

## -------------------------------------------------------------------------- ##

## load packages 
library(dplyr)
library(lubridate)
library(here)

## read in imputed tracks ## ------------------------------------------------- ##
imputes <- readRDS(here("pipeline","geoprocessed",
                        "all_behavlog_pseudotracks_20imputes_rerouted20mIso_geoprocessed_seafloor_2025Oct19.rds"))
imputes

# review
unique(imputes$DeployID)
summary(imputes)

# first: for analysis, we will consider dives as 50 meters or greater and 2 minutes
# or longer to maintain consistency across programming regimes. therefore, anything
# less than these thresholds will be considered a surface period 

# update the duration values based on timestamps 
beh <- imputes
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

## seafloor analysis prep ## ------------------------------------------------ ##

# filter dives + insular animals 
dives <- filter(beh2, What == "Dive") %>%
  filter(., DeployID != "PcTag090") %>%
  filter(., DeployID != "PcTag092") %>%
  filter(., DeployID != "PcTagP09") 

# assign a dive ID for each individual 
dives <- dives %>%
  group_by(DeployID, trk_id) %>%
  mutate(
    dive_id = row_number(),
    dive_id2 = paste0(DeployID,"-",dive_id)
  ) %>%
  ungroup()

# read in the dive IDs that should be included (4km error gap)
err_ids <- readRDS(here("pipeline","fkw_dive_ids_4kmerror_for_seafloor_analysis_2025Oct19.rds"))

# err_ids for insular animals (get tally)
err_ids %>%
  filter(., DeployID != "PcTag090") %>%
  filter(., DeployID != "PcTag092") %>%
  filter(., DeployID != "PcTagP09")  %>%
  tally()


# filter for those dives within the error criteria
dives_sub <- filter(dives, dive_id2 %in% err_ids$dive_id2)

# check that there's an equal number of dives per id and imputed track 
dives_chk <- dives_sub %>%
  group_by(DeployID, trk_id) %>%
  summarise(n = n())

# get the total number of dives per ID 
dives_id <- dives_sub %>%
  filter(., trk_id == 1) %>%
  group_by(DeployID) %>%
  tally()

# reassign object (easier)
dives <- dives_sub

# see if there are dives that are nans all around
nans <- filter(dives, is.na(depthH) & is.na(depthM) & is.na(depthL)) # nope

# create a new column for depth of insular animals, using the MHI and NWHI grids
# adjust depths of those just overlapping with -20m isobath grid cells. if point
# occurs in location where neither the MHI or NWHI grids had data, then use the
# GEBCO depth
dives <- dives %>%
  mutate(
    ins_depth = case_when(
      is.na(depthM) & !is.na(depthH) ~ depthH,
      is.na(depthM) & is.na(depthH) ~ depthL,
      is.na(depthH) & !is.na(depthM) ~ depthM,
      TRUE ~ depthH
    ),
    ins_slope = case_when(
      is.na(slopeM) & !is.na(slopeH) ~ slopeH,
      is.na(slopeM) & is.na(slopeH) ~ slopeL,
      is.na(slopeH) & !is.na(slopeM) ~ slopeM,
      TRUE ~ slopeH
    ),
    sfdepth = ifelse(ins_depth > -20, -20, ins_depth)
  )

summary(dives$ins_depth)
summary(dives$sfdepth)
summary(dives$ins_slope)

## isolate and summarize near seafloor dives ## ------------------------------ ##

# for each dive record, get summaries on seafloor depth across the 20 imputations
# and the differences between those summaries and the dive depth
dive_sf_sum <- dives %>%
  group_by(DeployID, dive_id2) %>%
  summarise(
    dive_depth = first(depth_avg50),
    dive_dur = first(dur_mins),
    dive_shape = first(Shape),
    tod = first(tod),
    mean_sf = mean(abs(sfdepth)),
    median_sf = median(abs(sfdepth)),
    sd_sf = sd(abs(sfdepth)),
    min_sf = min(abs(sfdepth)),
    max_sf = max(abs(sfdepth)),
    depth_diff_med = median_sf - dive_depth,
    depth_diff_mean = mean_sf - dive_depth,
    depth_diff_max = max_sf - dive_depth,
    prop_dive_sf = dive_depth/max_sf,
    mean_slope = mean(ins_slope),
    median_slope = median(ins_slope),
    sd_slope = sd(ins_slope),
    min_slope = min(ins_slope),
    max_slope = max(ins_slope)
  )

# save
write.csv(dive_sf_sum, here("outputs","summary_files",
                            "mhi_nwhi_dives_all_seafloor_summary_byID_bydive_2025Oct19.csv"), row.names = F)


# isolate dives that have standard deviation in seafloor depth < 100 m, and a
# maximum depth difference within 100m
dives_100sd <- filter(dive_sf_sum, sd_sf <= 100) %>%
  filter(., depth_diff_max <= 100)

dives_100sd_abs <- filter(dive_sf_sum, sd_sf <= 100) %>%
  filter(., abs(depth_diff_max) <= 100)

# four dives that are excluded when using absolute value. upon review, one of these
# dives fits criteria (others are more uncertain); include this one 
dives_100sd_abs <- bind_rows(dives_100sd_abs, dives_100sd[5,])

# save
write.csv(dives_100sd_abs, here("outputs","summary_files",
                            "mhi_nwhi_dives_sd100m_absdepthdiff100m_seafloor_summary_byID_bydive_2025Oct19.csv"), row.names = F)


# summarize the number of these dives by ID 
dd_sd100_sum <- dives_100sd_abs %>%
  group_by(DeployID) %>%
  summarise(
    n_sd100 = n(),
    mean_dd = mean(dive_depth),
    median_dd = median(dive_depth),
    sd_dd = sd(dive_depth),
    min_dd = min(dive_depth),
    max_dd = max(dive_depth),
    mean_prop = mean(prop_dive_sf),
    median_prop = median(prop_dive_sf),
    sd_prop = sd(prop_dive_sf),
    min_prop = min(prop_dive_sf),
    max_prop = max(prop_dive_sf),
    mean_du = mean(dive_dur),
    median_du = median(dive_dur),
    min_du = min(dive_dur),
    max_du = max(dive_dur),
    grmean_sf = mean(mean_sf),
    grmedian_sf = median(median_sf),
    grsd_sf = sd(sd_sf),
    grmin_sf = min(min_sf),
    grmax_sf = max(max_sf),
    mean_slope = mean(mean_slope),
    median_slope = median(median_slope),
    min_slope = min(min_slope),
    max_slope = max(max_slope)
  )
write.csv(dd_sd100_sum, here("outputs","summary_files",
                             "mhi_nwhi_probable_seafloor_dives_summary_byID_2025Oct19.csv"), row.names = F)

# get mean of medians 
mean(dd_sd100_sum$median_dd)
mean(dd_sd100_sum$median_du)

# get the number of dives within SD 100 by dive shape 
dd_sd100_shp <- dives_100sd_abs %>%
  group_by(DeployID, dive_shape) %>%
  summarise(
    n = n()
  ) %>%
  tidyr::pivot_wider(names_from = dive_shape, values_from = n) %>%
  dplyr::select(DeployID, Square, U, V)


# summary by time of day
dd_dive_sum_tod <- dives_100sd_abs %>%
  mutate(
    tod = factor(tod, levels = c("dawn","day","dusk","night"))
  ) %>%
  group_by(DeployID, tod) %>%
  summarise(
    n = n()
  ) %>%
  tidyr::pivot_wider(names_from = tod, values_from = n) %>%
  dplyr::select(DeployID, dawn, day, dusk, night)

# get the long format 
dd_dive_sum_tod_l <- dives_100sd_abs %>%
  mutate(
    tod = factor(tod, levels = c("dawn","day","dusk","night"))
  ) %>%
  group_by(DeployID, tod) %>%
  summarise(
    n_sf = n()
  )

# read in the total dive rate metrics by diel period (time of day/tod)
merge_tod <- readRDS(here("pipeline","dive_rates_by_tod_by_id_2025Oct10.rds"))

# estimate seafloor dive rate by diel period using the total hours of behavior
# data
merge_tod2 <- left_join(merge_tod, dd_dive_sum_tod_l, by = c("DeployID","tod")) %>%
  mutate(
    n_sf = ifelse(is.na(n_sf), 0, n_sf),
    sf_rate = n_sf/tot_hrs
  ) %>%
  filter(., DeployID != "PcTag090") %>%
  filter(., DeployID != "PcTag092") %>%
  filter(., DeployID != "PcTagP09")

# summarize average rates across individuals 
merge_tod_sum <- merge_tod2 %>%
  group_by(tod) %>%
  summarise(
    mean_sf_rate = mean(sf_rate),
    sd_sf_rate = sd(sf_rate),
    cv_sf_rate = sd_sf_rate/mean_sf_rate,
    min_sf_rate = min(sf_rate),
    max_sf_rate = max(sf_rate)
  ) 

# save objects for plotting later 
save(dive_sf_sum, dives_100sd_abs, merge_tod2, imputes, dives, file =
     here("pipeline","data_objects_for_seafloor_dives_analysis_2025Oct19.RData"))
