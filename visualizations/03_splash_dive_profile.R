## 03_splash_dive_profile.R: make dive profile plot (dive depth ~ time) of 
## SPLASH-tagged false killer whale for paper

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 31 Oct 2025

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

## source helper functions 
source(here("code","visualizations", "dive_plot_all_function_generic.R"))
source(here("code","visualizations", "dive_plot_day_function.R"))

## set tag ID 
Tag = "PcTag099"

## read in files for focal tag ## -------------------------------------------- ##

# behavior log
behavior <- read.csv("data/behavior/PcTag099_261474-Behavior.csv", header=T, stringsAsFactors = F)

# douglas filtered location data (need locs for getting sunrise/sunset times)
LocsFile <- read.csv("data/location/PcTag001-099_DouglasFiltered_r20d3lc2_ArgosGPS_2025OCTv1.csv", header = T)
locs <- filter(LocsFile, animal == "PcTag099") # filter out focal tag

# format datetimes: pay attention to the name of the date column in the Douglas filtered file 
locs$dateUTC <- as.POSIXct(locs$date, tz = "UTC")
locs$dateHST <- with_tz(locs$dateUTC, tzone = "Pacific/Honolulu")
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
sunriseset_sub <- sunriseset[c(1:33),]

# re-assign to sunriseset
sunriseset <- sunriseset_sub

## Specify xmin and xmax for "nights" in function (shaded boxes). 
## The first value/day for sunrise times needs to be removed for the math
## to work correctly. 

# sunset
xmin <- sunriseset$sunset
tz(xmin) # check time zone

# sunrise 
xmax <- sunriseset$sunrise
tz(xmax) # check time zone

# remove first xmax so shading math works correctly, so this will be the 
# second row (2) through the last row of the sunriseset data frame (16 in this
# example)
xmax_sub <- xmax[c(2:nrow(sunriseset))]

# add NA for last xmax, but still needs to be a POSIXct object
na_time <- with_tz(as.POSIXct(NA), tzone="Pacific/Honolulu")
tz(na_time) # check 

# bind back together with xmax
xmax_v2 <- append(xmax_sub, na_time)
xmax_v2 # check

# re-assign for simplicity
xmax <- xmax_v2

# create nights dataframe for plot function
nights = data.frame(xmin, xmax)

### plot dives over entire deployment ###

# specify folder to save plots to (set this to wherever you want plots to save) #
plot_folder = "outputs/dive_profile_plots/"

## plot for entire duration ## ----------------------------------------------- ##
## this requires: behavior file, title (optional), gaps (true or false to be included),
## plot color, depth limit (make sure is more than maximum dive depth in behavior log),
## ***double check file time zone (filetz)***, option for grid, color of gap blocks,
## depth axis labels and breaks. 
## Other options: see function arguments above

# get the max depth to inform the y axis limits
max(behavior$DepthMax, na.rm = T)

# make the plot
all <- DivePlotAll(behavior, title = "", gaps = TRUE, plotcolor = "#015b58", depthlim = -1200,
                   datebreaks = "1 day", nights = nights, filetz = "UTC", grid = FALSE,
                   plottz = "Pacific/Honolulu",
                   gapcolor = "black", depthlab = as.character(seq(-1200, 0, by = 100)),
                   depths = seq(-1200, 0, by = 100),
                   lineweight = 0.4)

all

# add aesthetic details 
all <-  all + theme(axis.text = element_text(size = 11, colour = "black"),
                    axis.title = element_text(size = 12, face = "bold"),
                    #axis.text.x = element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1)) + xlab("") + ylab("")

# plot it 
all

# save it 
ggsave(paste0(plot_folder, Tag, "_all_dives_2025Oct30.png"), all, device = "png",
       units = "in", width = 8, height = 4.75, dpi = 300)

## plot dives for a specific period ## --------------------------------------- ##
# change the datelim = to the two dates between 24 hour period
day1 <- DivePlotDay(behavior, title = "", gapcolor = "black", plotcolor = "#015b58", depthlim = -1200,
                    datebreaks = "6 hour", nights = nights,
                    datelim = as.POSIXct(c("2025-08-04 12:00:00", "2025-08-09 12:00:00"), tz = "Pacific/Honolulu"),
                    filetz = "UTC",
                    plottz = "Pacific/Honolulu",
                    depthlab = c(0,-250,-500,-750,-1000),
                    depths = c(0,-250,-500,-750,-1000),
                    gaps = F,
                    grid = FALSE,
                    lineweight = 0.4)

# change aesthetics
day1 <- day1 + theme(axis.text.y = element_text(size = 13, colour = "black"),
                     axis.title = element_text(size = 14),
                     axis.text = element_text(size = 12, colour = "black"),
                     #axis.text.x = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=1))

# change aesthetics
day1 <- day1 + xlab("Time (HST)") + ylab("Dive depth (m)")

# plot 
day1

# save it 
ggsave(paste0(plot_folder, Tag, "_Dives", "_2025Aug04-09.pdf"), day1,  units = "in", width = 10,
       height = 5, dpi = 300)
