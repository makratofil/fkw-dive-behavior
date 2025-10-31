## geoprocess.R: populate behavior log pseudotracks with temporal and spatial
## variables for further analysis.

## 03a: regular pseudotracks (i.e., single track, not multiple imputation)

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 09 Oct 2025

## ---------------------------------------------------------------------------- ##

## load packages 
library(lubridate)
library(oce)
library(lunar)
library(suntools)
library(sf)
library(stars)
library(purrr)
library(raster)
library(dplyr)
library(here)

## read in all behavior log pseudotracks and prep for geoprocessing ## -------- ##
files <- list.files(path = here("pipeline","crawl_pseudotracks"),
                    pattern = "2025Oct04.csv", full.names = T, recursive = F)
files

# function to read in and format the files 
fread <- function(x){
  
  d <- read.csv(x, header = T, stringsAsFactors = F)
  d$start_utc <- as.POSIXct(d$start_utc, tz = "UTC")
  d$start_hst <- as.POSIXct(d$start_hst, tz = "Pacific/Honolulu")
  
  return(d)
}

# apply the function
dfs <- lapply(files, fread)
tags <- bind_rows(dfs)
str(tags)

# subset the new tags (we already have some information for the old tags)
tags <- filter(tags, DeployID %in% c("PcTag095","PcTag096","PcTag097","PcTag099"))

## raster-based variables ## ------------------------------------------------- ##

# source the geoprocessing helper functions (specifically, devel version of 
## oce functions)
source(here("code","data_processing","oce_sun_moon_angle_devel_functions.R"))

# read in raster files (three different resolutions). read in from outside 
# directory so don't have to unnecessarily copy large files into this 
# directory

rasts <- "D:/04_OSU/06_Thesis/05_Analyses/pc_tag_processing/data/rasters/"

depthL <- stars::read_stars(paste0(rasts, "GebcodepthWUTM4N.nc"))

depthM <- stars::read_stars(paste0(rasts, "FalkorDepthUTM4N.nc"))

depthH <- stars::read_stars(paste0(rasts, "multibeamUTM4N.nc"))

slopeL <- stars::read_stars(paste0(rasts, "GebcoSlopeWUTM4N.nc"))

slopeM <- stars::read_stars(paste0(rasts, "FalkorSlopeUTM4N.nc"))

slopeH <- stars::read_stars(paste0(rasts, "MultibeamSlopeUTM4N.nc"))

# make tags dataframe into a spatial object 
tags_sf <- st_as_sf(tags, coords = c("lon","lat"), crs = 4326) %>%
  st_transform(crs = st_crs(depthH))

tags_sf <- tags_sf %>%
  mutate(
    # bathymetric variables
    depthH = as.numeric(stars::st_extract(depthH, ., method = "bilinear")[[1]]),
    slopeH = as.numeric(stars::st_extract(slopeH, ., method = "bilinear")[[1]]),
    depthM = as.numeric(stars::st_extract(depthM, ., method = "bilinear")[[1]]),
    slopeM = as.numeric(stars::st_extract(slopeM, ., method = "bilinear")[[1]]),
    depthL = as.numeric(stars::st_extract(depthL, ., method = "bilinear")[[1]]),
    slopeL = as.numeric(stars::st_extract(slopeL, ., method = "bilinear")[[1]])
  ) %>%
  # for temporal variables, get back into geographic coords
  st_transform(crs = 4326) %>%
  mutate(
    sunrise = as.POSIXct(suntools::sunriset(., start_hst, direction = "sunrise", POSIXct.out=T)$time),
    sunset = as.POSIXct(suntools::sunriset(., start_hst, direction = "sunset", POSIXct.out=T)$time),
    solar_noon = as.POSIXct(suntools::solarnoon(., start_hst, POSIXct.out=T)$time),
    civil_dawn = as.POSIXct(suntools::crepuscule(., start_hst, solarDep = 6, direction="dawn", POSIXct.out=T)$time),
    end_dawn = as.POSIXct(suntools::crepuscule(., start_hst, solarDep = -6, direction="dawn", POSIXct.out=T)$time),
    civil_dusk = as.POSIXct(suntools::crepuscule(., start_hst, solarDep = 6, direction="dusk", POSIXct.out=T)$time),
    start_dusk = as.POSIXct(suntools::crepuscule(., start_hst, solarDep = -6, direction="dusk", POSIXct.out=T)$time)
  )

summary(tags_sf)
str(tags_sf)


## add oceanographic variables ## -------------------------------------------- ##

# read in the old files 
old <- readRDS(here("pipeline","geoprocessed","all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_2025Feb18.rds"))

# make into sf object
old_sf <- st_as_sf(old, coords = c("lon","lat"), crs = 4326) %>%
  st_transform(crs = st_crs(tags_sf))

# bind with new tags
tags_sf <- bind_rows(old_sf, tags_sf)
class(tags_sf)

# create a deployment period variable for each tag that aligns with the variable
# files 
tags_sf <- tags_sf %>%
  mutate(
    dep_period = case_when(
      DeployID %in% c("PcTag026","PcTag028","PcTag030","PcTag032") ~ "2010",
      DeployID %in% c("PcTag035") ~ "2012",
      DeployID %in% c("PcTag037") ~ "2013",
      DeployID %in% c("PcTag049") ~ "2015",
      DeployID %in% c("PcTag055") ~ "2017",
      DeployID %in% c("PcTag074") ~ "2021",
      DeployID %in% c("PcTag090","PcTag092") ~ "2023",
      DeployID %in% c("PcTagP09") ~ "2024",
      DeployID %in% c("PcTag095","PcTag096","PcTag097") ~ "2025a",
      DeployID %in% c("PcTag099") ~ "2025b"
    )
  )

# check
unique(tags_sf$dep_period)

## create list of netCDF files to extract 

# daily chlorophyll-a
chd.files <- list.files(here("data","copernicus_oceancolor_daily")) 
chd.files

# monthly chlorophyll-a
chm.files <- list.files(here("data","copernicus_oceancolor_monthly")) 
chm.files

# all other oceanographic variables
oc.files <- list.files(here("data","copernicus_global_ocean_physics_reanalysis_all")) 
oc.files


# read in each daily chla file for each deployment year 
chlad2010 <- read_ncdf(here("data","copernicus_oceancolor_daily", chd.files[1]), eps = 1e-3)
chlad2012 <- read_ncdf(here("data","copernicus_oceancolor_daily", chd.files[2]), eps = 1e-3)
chlad2013 <- read_ncdf(here("data","copernicus_oceancolor_daily", chd.files[3]), eps = 1e-3)
chlad2015 <- read_ncdf(here("data","copernicus_oceancolor_daily", chd.files[4]), eps = 1e-3)
chlad2017 <- read_ncdf(here("data","copernicus_oceancolor_daily", chd.files[5]), eps = 1e-3)
chlad2021 <- read_ncdf(here("data","copernicus_oceancolor_daily", chd.files[6]), eps = 1e-3)
chlad2023 <- read_ncdf(here("data","copernicus_oceancolor_daily", chd.files[7]), eps = 1e-3)
chlad2024 <- read_ncdf(here("data","copernicus_oceancolor_daily", chd.files[8]), eps = 1e-3)
chlad2025a <- read_ncdf(here("data","copernicus_oceancolor_daily", chd.files[9]), eps = 1e-3)
chlad2025b <- read_ncdf(here("data","copernicus_oceancolor_daily", chd.files[10]), eps = 1e-3)

# read in each monthly chla file for each deployment year 
chlam2010 <- read_ncdf(here("data","copernicus_oceancolor_monthly", chm.files[1]), eps = 1e-3)
chlam2012 <- read_ncdf(here("data","copernicus_oceancolor_monthly", chm.files[2]), eps = 1e-3)
chlam2013 <- read_ncdf(here("data","copernicus_oceancolor_monthly", chm.files[3]), eps = 1e-3)
chlam2015 <- read_ncdf(here("data","copernicus_oceancolor_monthly", chm.files[4]), eps = 1e-3)
chlam2017 <- read_ncdf(here("data","copernicus_oceancolor_monthly", chm.files[5]), eps = 1e-3)
chlam2021 <- read_ncdf(here("data","copernicus_oceancolor_monthly", chm.files[6]), eps = 1e-3)
chlam2023 <- read_ncdf(here("data","copernicus_oceancolor_monthly", chm.files[7]), eps = 1e-3)
chlam2024 <- read_ncdf(here("data","copernicus_oceancolor_monthly", chm.files[8]), eps = 1e-3)
chlam2025a <- read_ncdf(here("data","copernicus_oceancolor_monthly", chm.files[9]), eps = 1e-3)
chlam2025b <- read_ncdf(here("data","copernicus_oceancolor_monthly", chm.files[10]), eps = 1e-3)

# read in each oc file for each deployment year
oc2010_3d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[1]), eps = 1e-3,
                       var = c("so","thetao","uo","vo")) %>%
  slice(., index = 1, along = "depth")
oc2010_2d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[1]), eps = 1e-3,
                       var = c("mlotst","zos"))

oc2012_3d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[2]), eps = 1e-3,
                    var = c("so","thetao","uo","vo")) %>%
  slice(., index = 1, along = "depth")
oc2012_2d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[2]), eps = 1e-3,
                       var = c("mlotst","zos"))

oc2013_3d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[3]), eps = 1e-3,
                       var = c("so","thetao","uo","vo")) %>%
  slice(., index = 1, along = "depth")
oc2013_2d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[3]), eps = 1e-3,
                       var = c("mlotst","zos"))

oc2015_3d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[4]), eps = 1e-3,
                       var = c("so","thetao","uo","vo")) %>%
  slice(., index = 1, along = "depth")
oc2015_2d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[4]), eps = 1e-3,
                       var = c("mlotst","zos"))

oc2017_3d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[5]), eps = 1e-3,
                       var = c("so","thetao","uo","vo")) %>%
  slice(., index = 1, along = "depth")
oc2017_2d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[5]), eps = 1e-3,
                       var = c("mlotst","zos"))

oc2021_3d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[6]), eps = 1e-3,
                       var = c("so","thetao","uo","vo")) %>%
  slice(., index = 1, along = "depth")
oc2021_2d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[6]), eps = 1e-3,
                       var = c("mlotst","zos"))

oc2023_3d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[7]), eps = 1e-3,
                       var = c("so","thetao","uo","vo")) %>%
  slice(., index = 1, along = "depth")
oc2023_2d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[7]), eps = 1e-3,
                       var = c("mlotst","zos"))

oc2024_3d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[8]), eps = 1e-3,
                       var = c("so","thetao","uo","vo")) %>%
  slice(., index = 1, along = "depth")
oc2024_2d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[8]), eps = 1e-3,
                       var = c("mlotst","zos"))

oc2025a_3d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[9]), eps = 1e-3,
                        var = c("so","thetao","uo","vo")) %>%
  slice(., index = 1, along = "depth")
oc2025a_2d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[9]), eps = 1e-3,
                       var = c("mlotst","zos"))

oc2025b_3d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[10]), eps = 1e-3,
                        var = c("so","thetao","uo","vo")) %>%
  slice(., index = 1, along = "depth")
oc2025b_2d <- read_ncdf(here("data","copernicus_global_ocean_physics_reanalysis_all", oc.files[10]), eps = 1e-3,
                        var = c("mlotst","zos"))

# double check crs are the same between datasets
st_crs(oc2024_3d) == st_crs(chlad2024)

# nest the data by deployment year, so we can extract data per file  
tags_sf$dmy <- as.character(date(tags_sf$start_utc))
dmy_nest <- tags_sf %>%
  st_transform(crs = st_crs(chlad2010)) %>% # projection the same as the files (same for chla and mld)
  group_by(dmy) %>%
  tidyr::nest()

# create a function for identifying corresponding ocean vars files and extracting values at each
# of the points by date 
get_vars <- function(data){
  
  # for testing
  #data <- dmy_nest[[2]][[3]]
  
  # get the dmy 
  dmy <- as.character(date(first(data$start_utc)))
  
  # get the dep period 
  dp <- data$dep_period[1]
  
  # assign the dataset based on the deployment year 
  if(dp == "2010"){
    chlad.file <- chlad2010
    chlam.file <- chlam2010
    oc3d.file <- oc2010_3d
    oc2d.file <- oc2010_2d
  } else if(dp == "2012"){
    chlad.file <- chlad2012
    chlam.file <- chlam2012
    oc3d.file <- oc2012_3d
    oc2d.file <- oc2012_2d
  } else if(dp == "2013"){
    chlad.file <- chlad2013
    chlam.file <- chlam2013
    oc3d.file <- oc2013_3d
    oc2d.file <- oc2013_2d
  } else if(dp == "2015"){
    chlad.file <- chlad2015
    chlam.file <- chlam2015
    oc3d.file <- oc2015_3d
    oc2d.file <- oc2015_2d
  } else if(dp == "2017"){
    chlad.file <- chlad2017
    chlam.file <- chlam2017
    oc3d.file <- oc2017_3d
    oc2d.file <- oc2017_2d
  } else if(dp == "2021"){
    chlad.file <- chlad2021
    chlam.file <- chlam2021
    oc3d.file <- oc2021_3d
    oc2d.file <- oc2021_2d
  } else if(dp == "2023"){
    chlad.file <- chlad2023
    chlam.file <- chlam2023
    oc3d.file <- oc2023_3d
    oc2d.file <- oc2023_2d
  } else if(dp == "2024"){
    chlad.file <- chlad2024
    chlam.file <- chlam2024
    oc3d.file <- oc2024_3d
    oc2d.file <- oc2024_2d
  } else if(dp == "2025a"){
    chlad.file <- chlad2025a
    chlam.file <- chlam2025a
    oc3d.file <- oc2025a_3d
    oc2d.file <- oc2025a_2d
  } else if(dp == "2025b"){
    chlad.file <- chlad2025b
    chlam.file <- chlam2025b
    oc3d.file <- oc2025b_3d
    oc2d.file <- oc2025b_2d
  }
  
  # filter the daily chla file for the dmy
  chlad_day <- filter(chlad.file, stringr::str_detect(time, as.character(dmy)))
  
  # get the date + 30 day lag and filter the raster for that 
  dmy30lag <- as.character(date(first(data$start_utc)) - days(30))
  chlad_30day <- filter(chlad.file, stringr::str_detect(time, as.character(dmy30lag)))
  
  # get the first day of the month and the previous month
  fm <- as.Date(format(as.Date(dmy), "%Y-%m-01"))
  fmp <- fm - days(30)
  
  # get monthly and 30-day lag monthly chla
  chlam_day <- filter(chlam.file, stringr::str_detect(time, as.character(fm)))
  chlam_30day <- filter(chlam.file, stringr::str_detect(time, as.character(fmp)))
  
  # filter the oc file for the dmy, for 3d variables 
  oc3d_day <- filter(oc3d.file, stringr::str_detect(time, as.character(dmy)))
  
  # filter the oc file for the dmy, for 2d variables
  oc2d_day <- filter(oc2d.file, stringr::str_detect(time, as.character(dmy)))
  
  # get spatial standard deviations of SST and SSH (across 3 pixels)
  sst_day <- oc3d_day["thetao"]
  ssh_day <- oc2d_day["zos"]
  
  sst_rast <- as(sst_day, "Raster")
  #plot(sst_rast)
  sst_agg <- aggregate(sst_rast, fact = 3, fun = sd)
  names(sst_agg) <- "sst_sd"
  #plot(sst_agg)
  sst_sd <- st_as_stars(sst_agg)
  
  ssh_rast <- as(ssh_day, "Raster")
  #plot(ssh_rast)
  ssh_agg <- aggregate(ssh_rast, fact = 3, fun = sd)
  names(ssh_agg) <- "ssh_sd"
  #plot(ssh_agg)
  ssh_sd <- st_as_stars(ssh_agg)
  
  # extract all values for each point within data
  data_new <- data %>%
    mutate(
      # daily contemporaneous chla
      chlad = as.numeric(stars::st_extract(chlad_day, .)$CHL),
      chlad_unc = as.numeric(st_extract(chlad_day, .)$CHL_uncertainty),
      chlad_flag = as.numeric(st_extract(chlad_day, .)$flags),
      # daily 30-day lag chla
      chla30d = as.numeric(stars::st_extract(chlad_30day, .)$CHL),
      chla30d_unc = as.numeric(st_extract(chlad_30day, .)$CHL_uncertainty),
      chla30d_flag = as.numeric(st_extract(chlad_30day, .)$flags),
      # monthly contemporaneous chla
      pp = as.numeric(stars::st_extract(chlam_day, .)$PP),
      pp_unc = as.numeric(st_extract(chlam_day, .)$PP_uncertainty),
      pp_flag = as.numeric(st_extract(chlam_day, .)$flags),
      # monthly 30-day lag chla
      # pp30 = as.numeric(stars::st_extract(chlam_30day, .)$PP),
      # pp30_unc = as.numeric(st_extract(chlam_30day, .)$PP_uncertainty),
      # pp30_flag = as.numeric(st_extract(chlam_30day, .)$flags),
      # all other dynamic vars
      sst = as.numeric(stars::st_extract(oc3d_day, .)$thetao),
      sst_sd = as.numeric(stars::st_extract(sst_sd, .)$sst_sd),
      ssh = as.numeric(stars::st_extract(oc2d_day, .)$zos),
      ssh_sd = as.numeric(stars::st_extract(ssh_sd, .)$ssh_sd),
      u_cur = as.numeric(stars::st_extract(oc3d_day, .)$uo),
      v_cur = as.numeric(stars::st_extract(oc3d_day, .)$vo),
      mld = as.numeric(stars::st_extract(oc2d_day, .)$mlotst)
    )
  
  # return the new dataframe
  return(data_new)
  
}

# apply the function
dmy_nest <- dmy_nest %>%
  mutate(
    oc_data = purrr::map(data, ~get_vars(.x))
  )

# unnest the data 
dmy_df <- dmy_nest %>%
  dplyr::select(dmy, oc_data) %>%
  tidyr::unnest() %>%
  sf::st_sf()

# get data back into regular dataframe with lat/lon, then add sun and moon
# angle and azimuth information (non-sf operation)
coords <- as.data.frame(st_coordinates(dmy_df))
pts_df <- bind_cols(dmy_df, coords) %>%
  rename(
    lon = X,
    lat = Y
  ) %>%
  arrange(DeployID, start_utc) %>%
  st_drop_geometry() %>%
  dplyr::select(-dmy) %>%
  mutate(
    ref = rep(T, length(start_utc)),
    sun_azimuth = sunAngle(start_utc, lon, lat, ref)$azimuth,
    sun_altitude = sunAngle(start_utc, lon, lat, ref)$altitude,
    moon_azimuth = moonAngle(start_utc, lon, lat, ref)$azimuth,
    moon_altitude = moonAngle(start_utc, lon, lat, ref)$altitude,
    moon_ill_fraction = moonAngle(start_utc, lon, lat, ref)$illuminatedFraction,
    moon_phase = lunar.phase(start_utc, name = T),
    tod = case_when(
      start_hst < civil_dawn | start_hst > civil_dusk ~ "night",
      start_hst > end_dawn & start_hst < start_dusk ~ "day",
      start_hst >= civil_dawn & start_hst <= end_dawn ~ "dawn",
      start_hst >= start_dusk & start_hst <= civil_dusk ~ "dusk"
    )
  )

# review the file 
summary(pts_df)

## save the file ## ---------------------------------------------------------- ##

# for geoprocessed to by split by time of day:
write.csv(pts_df, here("pipeline","geoprocessed","all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_2025Oct09.csv"), row.names = F)
saveRDS(pts_df, here("pipeline","geoprocessed","all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_2025Oct09.rds"))

