## 02_splash_data_summaries.R: data summaries for SPLASH tags

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 28 Apr 2025

## --------------------------------------------------------------------------- ##

## load packages
library(dplyr)
library(lubridate)
library(here)

## read in cleaned, processed behavior log data ## --------------------------- ##
beh <- readRDS(here("pipeline","clean_data_for_analysis",
                    "all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_split_tod_2025Feb18.rds"))

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

# add sex column 
beh <- beh %>%
  mutate(
    sex = case_when(
      DeployID %in% c("PcTag026","PcTag032","PcTag055","PcTag035","PcTag037","PcTagP09") ~ "Male",
      DeployID %in% c("PcTag090","PcTag092") ~ "Female",
      DeployID %in% c("PcTag028","PcTag030","PcTag074","PcTag049") ~ "Unknown"
    )
  )

# format datetime 
beh$start_utc <- as.POSIXct(beh$start_utc, tz = "UTC")
beh$end_utc <- as.POSIXct(beh$end_utc, tz = "UTC")

# read in fin length measurements and add to file 
fl <- read.csv(here("data","pc_splash_fin_length_for_R.csv"))

# join with data
beh <- left_join(beh,fl[,c(1,5)], by = "DeployID")

# make DeployID a factor leveled by population
beh$DeployID <- factor(beh$DeployID, levels = c("PcTag026","PcTag028","PcTag030","PcTag032","PcTag055",
                                                "PcTag074","PcTag035","PcTag037","PcTag049",
                                                "PcTag090","PcTag092","PcTagP09"))


## summary stats ## ---------------------------------------------------------- ##

# get the total number of dives that are 50m+ and 2min+
beh %>%
  group_by(What) %>%
  tally()

# get the number of days of behavior log data by tag ID
summary(beh$dur_days)
tot_days <- beh %>%
  group_by(DeployID) %>%
  summarise(
    tot_days = sum(dur_days),
    tot_hrs = sum(dur_hrs)
  )

sum(tot_days$tot_days) # 138.3
mean(tot_days$tot_days) # 11.5
median(tot_days$tot_days) # 9.7
range(tot_days$tot_days) # 3.8-20.1

# get the length of behavior log transmission duration (including gaps) for each
# ID, as to calculate data coverage
tot_length <- beh %>%
  group_by(DeployID) %>%
  summarise(
    trans_dur_sec = as.integer(as.duration(first(start_utc) %--% last(end_utc))),
    trans_dur_min = trans_dur_sec/60,
    trans_dur_hr = trans_dur_min/60,
    trans_dur_day = trans_dur_hr/24
  )

# get summary table of behavior log coverage by ID 
coverage <- left_join(tot_days, tot_length, by = "DeployID") %>%
  mutate(
    coverage = (tot_days/trans_dur_day)*100
  )


# get the total number of days by time of day category
tot_days_tod <- beh %>%
  group_by(tod) %>%
  summarise(
    tot_days = sum(dur_days)
  )

sum(tot_days_tod$tot_days) # check

# make a summary table by tag ID that shows the proportion of time at surface,
# the number of dives, etc. 
surf_days <- beh %>%
  filter(What == "Surface") %>%
  group_by(DeployID) %>%
  summarise(
    surf_days = sum(dur_days)
  )

n_dives <- beh %>%
  filter(What == "Dive") %>%
  group_by(DeployID) %>%
  summarise(
    n_dives = n()
  )

prop <- left_join(tot_days, surf_days, by = "DeployID") %>%
  left_join(., n_dives, by = "DeployID") %>%
  mutate(
    prop_surf = (surf_days/tot_days)*100,
    rate = n_dives/tot_hrs
  )

median(prop$rate) #0.39
mean(prop$rate) # 0.43
range(prop$rate) # 0.12 - 0.75
range(prop$prop_surf) # 89.5 - 99%


# make a summary table of dive summary stats by ID
dive_sum_id <- beh %>%
  filter(What == "Dive") %>%
  group_by(DeployID) %>%
  summarise(
    mean_depth = round(mean(depth_avg50),0),
    sd_depth = round(sd(depth_avg50),0),
    med_depth = round(median(depth_avg50),0),
    max_depth = round(max(depth_avg50),0),
    mean_dur = round(mean(dur_mins),1),
    sd_dur = round(sd(dur_mins),1),
    med_dur = round(median(dur_mins),1),
    max_dur = round(max(dur_mins),1)
  )

range(dive_sum_id$med_depth) # 112-744
range(dive_sum_id$med_dur) # 3.8-13.4


# get means of these metrics across median values for individuals 
dive_sum_means <- dive_sum_id %>%
  left_join(., prop[,c(1,7)], by = "DeployID") %>%
  summarise(
    mean_depth = mean(med_depth),
    cv_depth = sd(med_depth)/mean_depth,
    mean_dur = mean(med_dur),
    cv_dur = sd(med_dur)/mean_dur,
    mean_rate = mean(rate),
    cv_rate = sd(rate)/mean_rate
  )


# summarize by sex 
dive_sex_id_sum <- beh %>%
  filter(What == "Dive") %>%
  group_by(DeployID, sex) %>%
  summarise(
    mean_depth = mean(depth_avg50),
    sd_depth = sd(depth_avg50),
    med_depth = median(depth_avg50),
    max_depth = max(depth_avg50),
    mean_dur = mean(dur_mins),
    sd_dur = sd(dur_mins),
    med_dur = median(dur_mins),
    max_dur = max(dur_mins)
  ) %>%
  left_join(., prop[,c(1,7)], by = "DeployID") %>%
  group_by(sex) %>%
  summarise(
    mean_depth = mean(med_depth),
    cv_depth = sd(med_depth)/mean_depth,
    mean_dur = mean(med_dur),
    cv_dur = sd(med_dur)/mean_dur,
    mean_rate = mean(rate),
    cv_rate = sd(rate)/mean_rate
  )

# summarize by population
dive_pop_sum <- beh %>%
  filter(What == "Dive") %>%
  group_by(DeployID, population) %>%
  summarise(
    mean_depth = mean(depth_avg50),
    sd_depth = sd(depth_avg50),
    med_depth = median(depth_avg50),
    max_depth = max(depth_avg50),
    mean_dur = mean(dur_mins),
    sd_dur = sd(dur_mins),
    med_dur = median(dur_mins),
    max_dur = max(dur_mins)
  ) %>%
  left_join(., prop[,c(1,7)], by = "DeployID") %>%
  group_by(population) %>%
  summarise(
    mean_depth = mean(med_depth),
    cv_depth = sd(med_depth)/mean_depth,
    mean_dur = mean(med_dur),
    cv_dur = sd(med_dur)/mean_dur,
    mean_rate = mean(rate),
    cv_rate = sd(rate)/mean_rate
  )


# summarize within cluster 3
dive_c3_sum <- beh %>%
  filter(What == "Dive") %>%
  filter(DeployID %in% c("PcTag026","PcTag030","PcTag032","PcTag074")) %>%
  group_by(DeployID) %>%
  summarise(
    mean_depth = mean(depth_avg50),
    sd_depth = sd(depth_avg50),
    med_depth = median(depth_avg50),
    max_depth = max(depth_avg50),
    mean_dur = mean(dur_mins),
    sd_dur = sd(dur_mins),
    med_dur = median(dur_mins),
    max_dur = max(dur_mins)
  ) %>%
  left_join(., prop[,c(1,7)], by = "DeployID") %>%
  summarise(
    mean_depth = mean(med_depth),
    cv_depth = sd(med_depth)/mean_depth,
    mean_dur = mean(med_dur),
    cv_dur = sd(med_dur)/mean_dur,
    mean_rate = mean(rate),
    cv_rate = sd(rate)/mean_rate
  )


# subset dives and get overall values 
dives <- filter(beh, What == "Dive")
mean(dives$depth_avg50)
sd(dives$depth_avg50)
mean(dives$dur_mins)
sd(dives$dur_mins)

## summarize surface durations. first need to collapse consecutive durations - ## 
data <- beh
data$surface_id <- NA

# initialize variables
surface_id <- 0
last_dive <- FALSE

# iterate through the rows
for (i in 1:nrow(data)) {
  if (data$What[i] == 'Surface') {
    if (last_dive) {  # If last was a dive, we increment the surface_id
      surface_id <- surface_id + 1
      last_dive <- FALSE
    }
    data$surface_id[i] <- surface_id
  } else {
    last_dive <- TRUE  # When we encounter a Dive, reset last_dive
    data$surface_id[i] <- NA  # No surface_id for Dive rows
  }
}

# review
data_sub <- data %>%
  select(DeployID, What, surface_id, depth_avg50, dur_mins)

# summarize by indivdual 
surf_id_sum <- data %>%
  filter(What == "Surface") %>%
  group_by(DeployID, surface_id) %>%
  summarise(
    surf_dur = sum(dur_mins)
  ) %>%
  summarise(
    mean_dur = mean(surf_dur),
    median_dur = median(surf_dur),
    min_dur = min(surf_dur),
    max_dur = max(surf_dur),
    max_dur_days = max_dur/60/24,
    cv_dur = sd(surf_dur)/mean_dur
  ) 

range(surf_id_sum$median_dur)
mean(surf_id_sum$median_dur)

# summarize across individuals
surf_sum <- data %>%
  filter(What == "Surface") %>%
  group_by(DeployID, surface_id) %>%
  summarise(
    surf_dur = sum(dur_mins)
  ) %>%
  summarise(
    mean_dur = mean(surf_dur),
    median_dur = median(surf_dur),
    min_dur = min(surf_dur),
    max_dur = max(surf_dur),
    cv_dur = sd(surf_dur)/mean_dur
  ) %>%
summarise(
  gmean = mean(median_dur),
  cv = sd(median_dur)/gmean
)

# get the number of dives for spatial analyses
dives_space4km <- filter(beh, What == "Dive") %>%
  mutate(
    avg_se = (se.mu.x + se.mu.y)/2
  ) %>%
  filter(., avg_se < 4000) %>%
  group_by(DeployID) %>%
  tally()


## dive metrics by time of day summary stats ## ------------------------------- ##

# get dive rates and proportion of time at surface by time of day (table S3)
tot_dives_tod <- dives %>%
  group_by(DeployID, tod) %>%
  summarise(
    tot_dives = n()
  )

tot_surf_tod <- beh %>%
  filter(What == "Surface") %>%
  group_by(DeployID, tod) %>%
  summarise(
    tot_surf = sum(dur_hrs)
  )

tot_hr_tod <- beh %>%
  group_by(DeployID, tod) %>%
  summarise(
    tot_hrs = sum(dur_hrs)
  )


# get the proportion of dives by shape and time of day (table S3)
shape_tod <- dives %>%
  group_by(DeployID, tod, Shape) %>%
  summarise(
    n = n()
  ) %>%
  tidyr::pivot_wider(names_from = Shape, values_from = n) %>%
  mutate(
    Square = ifelse(is.na(Square),0,Square),
    U = ifelse(is.na(U), 0, U),
    V = ifelse(is.na(V), 0, V)
  )

# combine all
merge_tod <- left_join(tot_hr_tod, tot_dives_tod, by = c("DeployID", "tod")) %>%
  left_join(., tot_surf_tod, by = c("DeployID","tod")) %>%
  left_join(., shape_tod, by = c("DeployID","tod")) %>%
  mutate(tot_dives = ifelse(is.na(tot_dives), 0, tot_dives),
         rate = tot_dives/tot_hrs,
         prop_surf = (tot_surf/tot_hrs)*100,
         prop_u = (U/tot_dives)*100,
         prop_sq = (Square/tot_dives)*100,
         prop_v = (V/tot_dives)*100)

# save for later
saveRDS(merge_tod, here("pipeline","dive_rates_by_tod_by_id_2025Feb24.rds"))

# dive metrics by time of day (diel period)
dive_tod_sum <- beh %>%
  filter(What == "Dive") %>%
  group_by(DeployID, tod) %>%
  summarise(
    mean_depth = mean(depth_avg50),
    sd_depth = sd(depth_avg50),
    med_depth = median(depth_avg50),
    max_depth = max(depth_avg50),
    mean_dur = mean(dur_mins),
    sd_dur = sd(dur_mins),
    med_dur = median(dur_mins),
    max_dur = max(dur_mins)
  ) %>%
  left_join(., merge_tod, by = c("DeployID","tod")) %>%
  group_by(tod) %>%
  summarise(
    mean_depth = mean(med_depth),
    cv_depth = sd(med_depth)/mean_depth,
    mean_dur = mean(med_dur),
    cv_dur = sd(med_dur)/mean_dur,
    mean_rate = mean(rate),
    cv_rate = sd(rate)/mean_rate,
    mean_surf = mean(prop_surf),
    cv_surf = sd(prop_surf)/mean_surf,
    mean_sq = mean(prop_sq),
    cv_sq = sd(prop_sq)/mean_sq,
    mean_u = mean(prop_u),
    cv_u = sd(prop_u)/mean_u,
    mean_v = mean(prop_v),
    cv_v = sd(prop_v)/mean_v
  )

# get pearson's correlation between dive depth and duration
cor(dives$depth_avg50, dives$dur_mins, method = "spearman") # 0.72

## fin base length summaries ## ---------------------------------------------- ##

# pearson's correlation between median dive depth and fin base length
cor(dive_sum_id$med_depth, fl$fin_base, method = "spearman") # 0.36
cor(dive_sum_id$med_dur, fl$fin_base, method = "spearman") # 0.37

## save data objects for use in visualizations later ## ---------------------- ##
save(dives, merge_tod, file = here("pipeline","splash_data_summary_objects_for_figures_2025Apr29.RData"))

