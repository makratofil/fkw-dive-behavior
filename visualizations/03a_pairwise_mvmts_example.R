## 03_pairwise_mvmts_example.R: example/snapshot of pairwise movements (vertical
## and horizontal) on PcTag030-032 for manuscript figure 3

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 29 Apr 2025

## --------------------------------------------------------------------------- ##

## load libraries 
library(ggplot2)
library(lubridate)
library(gridExtra)
library(scales)
library(suncalc)
library(tidyr)
library(dplyr)
library(here)
library(ggpubr)
library(sf)
library(sfdep)

## load functions used to make the plots 
source(here("code","visualizations","dive_plot_day_function.R"))

## prep behavior log data for plotting ## ------------------------------------ ##

## read in qa/qc'd behaivor log file 
behavior <- readRDS(here("pipeline","qaqcd_behavior_data",
                         "PcTags_SPLASH_thruP09_behavior_logs_qaqcd_2025Feb17.rds"))
str(behavior) # quick review

# subset the individuals in the pair of interest
pair <- filter(behavior, DeployID %in% c("PcTag030","PcTag032"))

# standardize dives to 50m and 2mins
beh <- pair %>%
  mutate(DepthAvg50 = replace (DepthAvg, DepthAvg < 50, NA)) %>%
  mutate(What = ifelse(is.na(DepthAvg50), "Surface", "Dive"),
         dur_mins = DurationAvg/60) %>%
  mutate(depth_avg50 = ifelse (DepthAvg50 > 50 & dur_mins < 1.9, NA, DepthAvg50)) %>%
  mutate(What = 
           ifelse (is.na(DepthAvg50), "Surface", "Dive")) 

# assign the depthavg and durationavg columns to the original depth/duration min/max
# columns to align with how the plotting code works (it will all work out to plot
# the same)
pair <- beh
pair$DurationMin <- pair$DurationAvg
pair$DurationMax <- pair$DurationAvg
pair$DepthMin <- pair$DepthAvg50
pair$DepthMax <- pair$DepthAvg50

## read in douglas filtered data to get locations for sunrise/sunset times 
locations <- read.csv(here("data","location","PcTag001-092_DouglasFiltered_r20d3lc2_ArgosGPS_2024JUNv1.csv"), header = T)
locations$dateUTC <- as.POSIXct(locations$date, tz = "UTC")
locations$dateHST <- with_tz(locations$dateUTC, tzone = "Pacific/Honolulu")

# subset one of the pair (only need one, sunrise/sunset times will be the same 
# since they were tagged at the same time and moving together)
locs <- filter(locations, animal == "PcTag030") 

# format datetimes: pay attention to the name of the date column in the Douglas filtered file 
locs$date_ymd <- as.Date(locs$dateUTC, format = "%Y%m%d")
str(locs)
tz(locs$dateHST)

# compute average lat/lon for each day for suncalc
sunriseset <- locs %>%
  select(latitude, longitud, date_ymd) %>%
  group_by(date_ymd) %>%
  summarise(mean(latitude), mean(longitud)) %>%
  rename(lat = "mean(latitude)", lon = "mean(longitud)")

# change column name
colnames(sunriseset)[colnames(sunriseset) == "date_ymd"] <- "date"

# get sunrise/sunset times
sunriseset <- getSunlightTimes(data = sunriseset, keep = c("sunrise", "sunset"), tz = "Pacific/Honolulu")

# subset days where behavior data transmitted. for example, a tag may have transmitted
# LOCATION data from 8/8/2021 through 8/23/2021 but only BEHAVIOR data through 
# 8/18/2021, so subset the sunriseset data so we only plot the days with behavior
# data

# behavior start and end dates by ID 
behav_dates <- pair %>%
  group_by(DeployID) %>%
  summarise(
    start = first(start_utc),
    end = last(end_utc)
  )

# choose the latest start date (when both tags transmitting)
max(behav_dates$start)
behav_st <- as.Date(filter(behav_dates, start == max(start))$start)

# choose the earliest end date (when both tags transmitting)
behav_en <- as.Date(filter(behav_dates, end == min(end))$end)

# now filter the sunrise/sunset dataframe to only include these dates 
sunriseset_sub <- filter(sunriseset, between(date, behav_st, behav_en))

# Specify xmin and xmax for "nights" in function (shaded boxes). 
# The first value/day for sunrise times needs to be removed for the math
# to work correctly. 
# sunset
xmin <- sunriseset_sub$sunset
tz(xmin) # check time zone

# sunrise 
xmax <- sunriseset_sub$sunrise
tz(xmax) # check time zone

# remove first xmax so shading math works correctly, so this will be the 
# second row (2) through the last row of the sunriseset data frame 
xmax_sub <- xmax[c(2:nrow(sunriseset_sub))]

# add NA for last xmax, but still needs to be a POSIXct object
na_time <- with_tz(as.POSIXct(NA), tzone="Pacific/Honolulu")
tz(na_time) # check 

# bind back together with xmax
xmax_v2 <- append(xmax_sub, na_time)
xmax_v2 # check

# re-assign for simplicity
xmax <- xmax_v2

# create nights dataframe for plot function
nights <- data.frame(xmin, xmax)

## for mapping ## ------------------------------------------------------------ ##

# get CTCRW model data with calculated pairwise distances for mapping the tracks and 
# plotting proximity 
load(here("pipeline","crawl_pseudotracks","PcTag030_PcTag032_pairwise_dists_locs_2025Mar24v2.RData"))

# subset ids  
locs_id1_sf <- pair1 %>% st_transform(crs = 4326)

locs_id2_sf <- pair2 %>% st_transform(crs = 4326)

# get the pseudotracks for mapping points of dives 
pseu1 <- read.csv(here("pipeline","crawl_pseudotracks","PcTag030_behavlog_pseudotrack_rerouted20mIso_2025Feb18.csv")) %>%
  st_as_sf(., coords = c("lon","lat"), crs = 4326)
pseu1$start_utc <- as.POSIXct(pseu1$start_utc, tz = "UTC")
pseu1$end_utc <- as.POSIXct(pseu1$end_utc, tz = "UTC")

pseu2 <- read.csv(here("pipeline","crawl_pseudotracks","PcTag032_behavlog_pseudotrack_rerouted20mIso_2025Feb18.csv")) %>%
  st_as_sf(., coords = c("lon","lat"), crs = 4326)
pseu2$start_utc <- as.POSIXct(pseu2$start_utc, tz = "UTC")
pseu2$end_utc <- as.POSIXct(pseu2$end_utc, tz = "UTC")

# assign new dive/surface columns to pseudotracks 
id1_bl <- filter(pair, DeployID == "PcTag030")
pseu1$What <- id1_bl$What
id2_bl <- filter(pair, DeployID == "PcTag032")
pseu2$What <- id2_bl$What

# read in coastline for mapping 
coast <- st_read(here("data","shapefiles","FisheriesIslands.shp"))

# ESRI ocean basemap
esri_ocean <- paste0('https://services.arcgisonline.com/arcgis/rest/services/',
                     'Ocean/World_Ocean_Base/MapServer/tile/${z}/${y}/${x}.jpeg')

## make the plots ## --------------------------------------------------------- ##

# specify folder to save plots to (set this to wherever you want plots to save) #
plot_folder = here("outputs","dive_profile_plots_pairs_2025Mar24","PcTag030_PcTag032")

# get the vector of dates to make plots for 
dates <- sunriseset_sub$date
dates

# select the date of interest for this snapshot/example plot
d <- 22

# identify the span of dates to plot 
d1 <- dates[d-1]
d2 <- dates[d]
  
# get the HST times to subset (these are specific to the snapshot)
d1_plot <- paste0(d2, " 03:00:00")
d2_plot <- paste0(d2, " 19:00:00")
  
# identify the two individuals to plot 
id1 <- unique(pair$DeployID)[1]
id2 <- unique(pair$DeployID)[2]
  
## dive profile plots ## --------------------------------------------------- ##

# get the UTC dates to plot (data is in UTC, but we want to plot in HST)
d1_plot_utc <- as.POSIXct(d1_plot, tz = "UTC") + hours(10)
d2_plot_utc <- as.POSIXct(d2_plot, tz = "UTC") + hours(10)
pair_date <- filter(pair, start_utc >= d1_plot_utc & end_utc <= d2_plot_utc) 

# get the maximum dive depth out of both individuals for these dates for 
# setting the y-axis. make the y-axis upper limit 100m more than the max
max_depth <- -(round(max(pair_date$DepthMax, na.rm = T),-2) + 50)
  
  # make the plot for individual 1 
  id1_behav <- filter(pair, DeployID == id1)
  id1_plot <- DivePlotDay(id1_behav, title = "", gapcolor = "black", plotcolor = "#41476b",
                          depthlim = max_depth,
                          depths = seq(0,max_depth, by = -150),
                          as.character(c(seq(0,(max_depth), by = -150))),
                          datebreaks = "4 hour", nights = nights,
                          datelim = as.POSIXct(c(d1_plot, d2_plot), tz = "Pacific/Honolulu"),
                          filetz = "UTC",
                          plottz = "Pacific/Honolulu",
                          gaps = TRUE,
                          grid = FALSE, lineweight = 0.5) + 
    theme(
            #axis.title = element_text(size = 13),
            axis.text = element_text(size = 14, colour = "black"),
            axis.text.x = element_blank(),
            axis.title = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) 
  id1_plot
  ggsave(here(plot_folder,"PcTag030_PcTag032_id_profile_2011-01-02_for_ms_v2.png"), width = 7, height = 4, units = "in")
  
  # make the plot for individual 2 
  id2_behav <- filter(pair, DeployID == id2)
  id2_plot <- DivePlotDay(id2_behav, title = "", gapcolor = "black", plotcolor = "#c67b6f",
                          depthlim = max_depth,
                          depths = seq(0,max_depth, by = -150),
                          as.character(c(seq(0,(max_depth), by = -150))),
                          datebreaks = "4 hour", nights = nights,
                          datelim = as.POSIXct(c(d1_plot, d2_plot), tz = "Pacific/Honolulu"),
                          filetz = "UTC",
                          plottz = "Pacific/Honolulu",
                          gaps = TRUE,
                          grid = FALSE, lineweight = 0.5) + 
    theme(
      #axis.title = element_text(size = 13),
      axis.text = element_text(size = 14, colour = "black"),
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) 
  id2_plot
  ggsave(here(plot_folder,"PcTag030_PcTag032_id2_profile_2011-01-02_for_ms_v3.png"), width = 7, height = 4, units = "in")
  
## map of locations ## -------------------------------------------------------- ##
  
  # subset dive locations between the dates (locations also in UTC)
  id1_dives <- filter(pseu1, start_utc >= d1_plot_utc & end_utc <= d2_plot_utc) %>%
    filter(., What == "Dive") %>%
    mutate(timestamp = start_utc)
  id2_dives <- filter(pseu2, start_utc >= d1_plot_utc & end_utc <= d2_plot_utc) %>%
    filter(., What == "Dive") %>%
    mutate(timestamp = start_utc)
  
  # subset crawl locations between the dates, add pseudotracks, then make linestrings 
  id1_d_locs <- filter(locs_id1_sf, timestamp >= d1_plot_utc & timestamp <= d2_plot_utc) 
  
  id1_lines <- bind_rows(id1_d_locs, id1_dives) %>%
    arrange(timestamp) %>%
    summarise(do_union = F) %>%
    st_cast("LINESTRING")
  
  id2_d_locs <- filter(locs_id2_sf, timestamp >= d1_plot_utc & timestamp <= d2_plot_utc) 
  
  id2_lines <- bind_rows(id2_d_locs,id2_dives) %>%
    arrange(timestamp) %>%
    summarise(do_union = F) %>%
    st_cast("LINESTRING")
  
  # ellipses around locations 

    id1_ellipses <- id1_d_locs %>%
      st_transform(., crs = 26904) %>% # transform to UTM b/c error is in meters
      rowwise() %>%
      mutate(
        ellipse = st_ellipse(geometry = geometry, sx = se.mu.x, sy = se.mu.y)
      )
    
    id1_el <- id1_ellipses$ellipse
    id1_elps <- st_cast(id1_el, "POLYGON")

    
    id2_ellipses <- id2_d_locs %>%
      st_transform(., crs = 26904) %>% # tranform to UTM b/c error is in meters
      rowwise() %>%
      mutate(
        ellipse = st_ellipse(geometry = geometry, sx = se.mu.x, sy = se.mu.y)
      )
    
    id2_el <- id2_ellipses$ellipse
    id2_elps <- st_cast(id2_el, "POLYGON")

  # get bounding box for both
  s <- data.frame(id = c("id1","id2"),
                  xmin = c(st_bbox(id1_lines)[1], st_bbox(id2_lines)[1]),
                  xmax = c(st_bbox(id1_lines)[3], st_bbox(id2_lines)[3]),
                  ymin = c(st_bbox(id1_lines)[2], st_bbox(id2_lines)[2]),
                  ymax = c(st_bbox(id1_lines)[4], st_bbox(id2_lines)[4]))
  
  # set limits 
  xmin <- round(min(s$xmin - 0.05, na.rm = T),2)
  xmax <- round(max(s$xmax + 0.05, na.rm = T),2)
  ymin <- round(min(s$ymin - 0.2, na.rm = T),2)
  ymax <- round(max(s$ymax + 0.2, na.rm = T),2)
  

    # make the map
    map <- ggplot() +
      
      # background and coastline
      ggspatial::annotation_map_tile(type = esri_ocean, zoomin = 1, progress = "none") +
      geom_sf(data = coast, color = NA, fill = "gray19") +
      
      # ellipses
      geom_sf(data = id1_elps, fill = "#41476b", color = "#41476b", alpha = 0.2, linetype = "dashed") +
      geom_sf(data = id2_elps, fill = "#c67b6f", color = "#c67b6f", alpha = 0.2, linetype = "dashed") +

      # tracklines
      geom_sf(data = id1_lines, color = "#41476b") +
      geom_sf(data = id2_lines, color = "#c67b6f") +

      # dives
      geom_sf(data = id1_dives, shape = 21, color = "#41476b", fill = "white", size = 1.75) +
      geom_sf(data = id2_dives, shape = 21, color = "#c67b6f", fill = "white", size = 1.75) +
      
      # start locations of the timeseries for both tags
      geom_sf(data = id1_d_locs[1,], shape = 24, color = "#41476b",fill = "white", size = 3) +
      geom_sf(data = id2_d_locs[1,], shape = 24, color = "#c67b6f",fill = "white", size = 3) +
      
      # aesthetics
      coord_sf(crs = 4326,
               xlim = c(xmin, xmax),
               ylim = c(ymin, ymax)) +
      scale_y_continuous(breaks = c(20.2,20.4,20.6,20.8)) +
      theme_bw() +
      theme(
        axis.text = element_text(color = "black", size = 8),
        axis.ticks.length=unit(-0.1, "cm"),
        legend.position = "none"
      ) +
      ggspatial::annotation_scale(text_cex = 0.8, text_col = "black") 
    
    map
    
    ggsave(here(plot_folder,"PcTag030_PcTag032_pair_map_2011-01-02_for_ms_v7.png"), width = 7, height = 5, units = "in")

## distance between pairs plot ## -------------------------------------------- ##
    
    # subset timeframe
    dist <- filter(pair_final, timestamp >= d1_plot_utc & timestamp <= d2_plot_utc)
    
    # add dates in HST
    dist$timestamp_hst <- as.POSIXct(format(dist$timestamp, tz="Pacific/Honolulu"), tz="Pacific/Honolulu")
    id1_dives$start_hst <- as.POSIXct(id1_dives$start_hst, tz = "Pacific/Honolulu")
    id2_dives$start_hst <- as.POSIXct(id2_dives$start_hst, tz = "Pacific/Honolulu")
    
    d1_plot_hst <- as.POSIXct(format(d1_plot_utc, tz = "Pacific/Honolulu"), tz = "Pacific/Honolulu")
    d2_plot_hst <- as.POSIXct(format(d2_plot_utc, tz = "Pacific/Honolulu"), tz = "Pacific/Honolulu")

    
    ggplot(dist, aes(x = timestamp_hst, y = dist_km)) +
      
      # shaded rectangles for nighttime periods
      annotate(geom = "rect", xmin = d1_plot_hst, xmax = nights$xmax[21], ymin = -Inf, ymax = Inf,
               fill = "grey", alpha = 0.5, color = NA) +
      annotate(geom = "rect", xmin = nights$xmin[22], xmax = d2_plot_hst, ymin = -Inf, ymax = Inf,
               fill = "grey", alpha = 0.5, color = NA) +
      
      # horizontal line for distance b/t pair = 0
      geom_line(aes(y = 0), linetype = "dashed") +
      
      # line for distance b/t pair over time
      geom_line() + 
      
      # horizontal lines + points for plotting the timeseries of dives below 0
      geom_line(aes(y = -0.25), color = "#41476b") +
      geom_line(aes(y = -0.5), color = "#c67b6f") +
      geom_point(data = id1_dives, aes(x = start_hst, y = -.25), color = "#41476b") +
      geom_point(data = id2_dives, aes(x = start_hst, y = -.5), color = "#c67b6f") +
 
      # aesthetics
      scale_x_datetime(limits = c(d1_plot_hst,d2_plot_hst),
                       labels = date_format("%H:%M", tz = "Pacific/Honolulu"),
                       expand = c(0,0),
                       breaks = "4 hours") +
      scale_y_continuous(breaks = c(0,1,2,3,4,5)) +
      ylab("") +
      theme_bw() +
      theme(
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 14, color = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
      ) 

    ggsave(here(plot_folder,"PcTag030_PcTag032_pair_dist_2011-01-02_for_ms_v8.png"), width = 7, height = 4, units = "in")
