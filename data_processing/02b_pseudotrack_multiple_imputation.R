## pseudotrack_multiple_imputation.R: use fitted crawl models to simulate 
## multiple imputations of the behavior log "pseudotrack" conditioned on 
## model MLE parameter estimates

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 23 Feb 2025

## --------------------------------------------------------------------------- ##

## load package libraries
library(dplyr)
library(here)
library(lubridate)
library(readr)
library(purrr)
library(sf)
library(crawl) 
library(ggplot2)
library(janitor)
library(zoo)
library(pathroutr)

## read in the fitted movement models from "pseudotrack.R" ## ---------------- ##
mov <- readRDS(here("pipeline","PcTags_SPLASH_thru092_Kalman_Crawl_Fitted_Model_wSegments_2025Feb18.rds"))

## get objects for pathroutr to reroute predicted locations around land ------ ##

# 20m isobath and bathymetry for augmented points
iso20 <- readRDS(here("code","data_processing","isobath20m_for_rerouting.rds"))
bathy <- readRDS(bathy, here("code","data_processing","raster_for_rerouting.rds"))

# quick map
ggplot()+
  geom_sf(data=iso20)+
  theme_classic()+ theme(legend.position="none")

# check crs
st_crs(iso20)

# Visability graph for rerouting paths around barrier (iso20)
iso20_b <- st_buffer(iso20, 10000)

# add augmented points for isobaths that aren't as complete, like niihau, Kauai,
# nihoa, and necker
aug_pts <- raster::rasterToContour(bathy, levels=c(-100, -200, -300)) %>% st_as_sf() %>% 
  st_sample(300, type='regular') %>% st_cast("POINT")
aug_pts <- aug_pts[st_intersects(aug_pts, iso20_b, sparse=F)]  

# quick check of aug points 
ggplot() +
  geom_sf(data = iso20) +
  geom_sf(data = aug_pts)

# create visibility graph
vis_graph <- pathroutr::prt_visgraph(iso20, centroids=TRUE, aug_points=aug_pts)
#plot(vis_graph)

# going to loop through each animal, read in the behavior log file, and provide
# the timestamps of the behavior log file to crawl to predict locations at 
# (this is a simplified approach to the pseudotracks)
ids <- unique(mov$animal)

for(i in 1:length(ids)){
  # for testing
  #i = 12
  
  # select the animal
  id <- ids[i]
  
  # get the compiled QA/QC'd behavior log file 
  bq <- readRDS(here("pipeline","qaqcd_behavior_data",
                     "PcTags_SPLASH_thruP09_behavior_logs_qaqcd_2025Feb17.rds"))
  
  # filter the behavior log for the ID  
  id_dive <- filter(bq, DeployID == id)
  
  # format the datetime columns according to time zone it was programmed in,
  # and make a standardized datetime_utc column 
  if(id %in% c("PcTag035","PcTag037")){
    TZ = "Pacific/Honolulu"
  } else {
    TZ = "UTC"
  }
  
  id_dive <- id_dive %>%
    mutate(
      Start = as.POSIXct(gsub("\\.5", "", id_dive$Start), format = "%H:%M:%S %d-%b-%Y", tz = TZ),
      End = as.POSIXct(gsub("\\.5", "", id_dive$End), format = "%H:%M:%S %d-%b-%Y", tz = TZ),
      start_hst = as.POSIXct(format(Start, tz="Pacific/Honolulu"), tz="Pacific/Honolulu")
    )
  
  # create vector of times to predict crawl for
  log_times <- id_dive$start_utc
  
  # get the crawl fitted object for that ID
  id_crawl <- mov[[7]][[i]]
  
  # create simulation object from fitted model, using behavior log times as
  # the pred time
  tag_sim <- crawl::crwSimulator(id_crawl, parIS = 0, predTime = log_times)
  
  # reroute 20 imputations from simulated object above
  imputes <- tibble::tibble(.rows = 20) %>% 
    dplyr::rowwise() %>% 
    mutate(postis = list(crwPostIS(tag_sim, fullPost = FALSE)),
           pts = list(crw_as_sf(postis, locType = c("p","o"), ftype = "POINT")),
           pts = list(st_transform(pts, crs = 3857)),
           pts_trim = list(pathroutr::prt_trim(pts, iso20)),
           rrt_pts = list(pathroutr::prt_reroute(pts_trim, iso20, vis_graph, blend = FALSE)),
           pts_fix = list(pathroutr::prt_update_points(rrt_pts = rrt_pts, trkpts = pts_trim)),
           lines_fix = list(pts_fix %>% summarise(do_union = FALSE) %>% st_cast('LINESTRING'))
    )
  
  # get the number of points per simulated track 
  im1 <- imputes[[3]][[5]]
  n_pts <- nrow(im1)
  
  # get the simulated points in a dataframe, add a datetime_utc object, and then
  # add a track id for each imputation
  sim_points <- do.call(rbind, imputes$pts_fix)
  sim_points$start_utc <- sim_points$timestamp
  sim_points$trk_id = rep(1:ceiling(nrow(sim_points)/n_pts), each = n_pts)[1:nrow(sim_points)]
  
  # append the other variables from dive behavior log to the crawl file 
  sim_points2 <- left_join(sim_points, id_dive, by = "start_utc")

  # get the lat/lon coordinates back
  sim_points2_wgs <- st_transform(sim_points2, crs = 4326)
  id_coords <- as.data.frame(st_coordinates(sim_points2_wgs))
  
  # format/clean up the data 
  clean_sim <- sim_points2 %>%
    bind_cols(., id_coords) %>%
    st_drop_geometry() %>%
    dplyr::rename(
      lon = X,
      lat = Y
    ) %>%
    dplyr::select(DeployID, Ptt, timestamp, trk_id, locType, Start, End, What, Shape, DepthMin, DepthMax,
                  DurationMin, DurationMax, lon, lat,
                  start_utc, end_utc, start_hst) %>%
    rename(
      crawl_timestamp = timestamp
    )
  
  # save the file 
  write.csv(clean_sim, here("pipeline","crawl_pseudotracks_multiple_imputation",
                             paste0(id,"_behavlog_pseudotrack_20imputes_rerouted20mIso_2025Feb23.csv")), row.names = F)
  
}
