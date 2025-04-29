## pairwise_distances.R: calculate straight line distances between pairs of SPLASH
## tagged false killer whales 

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 22 Apr 2025

## --------------------------------------------------------------------------- ##

## load packages 
library(dplyr)
library(here)
library(lubridate)
library(sf)
library(crawl)

## note: this script is set up to run one pair at a time through the workflow

## read in behavior log data of pairs ## ------------------------------------- ##
beh <- readRDS(here("pipeline","clean_data_for_analysis",
                    "all_behavlog_pseudotracks_rerouted20mIso_geoprocessed_split_tod_2025Feb18.rds"))

# filter pair
pair <- filter(beh, DeployID %in% c("PcTag030","PcTag032"))
#pair <- filter(beh, DeployID %in% c("PcTag090","PcTag092"))

# get the start and end times by id
str(pair)
stends <- pair %>%
  group_by(DeployID) %>%
  summarise(
    start = first(start_utc),
    end = last(end_utc)
  )

# select the latest start time and the earliest end time to get start/end of 
# period of overlap
start <- filter(stends, start == max(start))$start
end <- filter(stends, end == min(end))$end

# get a sequence of times to estimate positions at during the period of behavior
# log overlap
times_split <- data.frame(start_utc = timeplyr::time_seq(start, end, time_by = "15 minutes"))

# get dive times too 
pair_dives <- data.frame(start_utc = filter(pair, What == "Dive")$start_utc)
pair_dives$what <- "dive"

# bind together 
times_all <- bind_rows(times_split, pair_dives) %>%
  arrange(start_utc)

## read in crawl fitted models and estimate positions at the specified times - ##
mov <- readRDS(here("pipeline","PcTags_SPLASH_thru092_Kalman_Crawl_Fitted_Model_wSegments_2024Jul01.rds"))
mov

# subset to the pair of interest
mov_sub <- filter(mov, animal %in% c("PcTag030","PcTag032"))

## get objects for pathroutr to reroute predicted locations around land 
# read in 20m isobath 
iso20 <- readRDS(here("code","data_processing","isobath20m_for_rerouting.rds"))
bathy <- readRDS(bathy, here("code","data_processing","raster_for_rerouting.rds"))

# Visability graph for rerouting paths around barrier (iso20)
iso20 <- st_transform(iso20, crs = 26904)
iso20_b <- st_buffer(iso20, 10000) %>% st_transform(crs = 26904)

# add augmented points for isobaths that aren't as complete, like niihau, Kauai,
# nihoa, and necker
aug_pts <- raster::rasterToContour(bathy, levels=c(-100, -200, -300)) %>% st_as_sf() %>% 
  st_sample(300, type='regular') %>% st_cast("POINT") %>% st_transform(crs = 26904)
aug_pts <- aug_pts[st_intersects(aug_pts, iso20_b, sparse=F)]  

# quick check of aug points 
ggplot() +
  geom_sf(data = iso20) +
  geom_sf(data = aug_pts)

# create visibility graph
vis_graph <- pathroutr::prt_visgraph(iso20, centroids=TRUE, aug_points=aug_pts)

# going to loop through each animal, read in the behavior log file, and provide
# the timestamps of the behavior log file to crawl to predict locations at 
# (this is a simplified approach to the pseudotracks)
ids <- unique(mov_sub$animal)

pseus <- list()

for(i in 1:length(ids)){
  # for testing
  #i = 1
  
  # select the animal
  id <- ids[i]
  
  # create vector of times to predict crawl for
  log_times <- times_all$start_utc
  
  # get the crawl fitted object for that ID
  id_crawl <- mov_sub[[7]][[i]]
  
  # predict crawl locations for the sequence of times
  id_pseu <- crwPredict(id_crawl, predTime = log_times, return.type = "flat") %>%
    crw_as_sf(.,  ftype = "POINT", locType = c("p","o")) %>%
    st_transform(., crs = 26904)
  
  # reroute the locations around land 
  pts_trim <- pathroutr::prt_trim(id_pseu, iso20)
  pts_rr <- pathroutr::prt_reroute(pts_trim, iso20, vis_graph, blend = F)
  pts_new <- pathroutr::prt_update_points(pts_rr, pts_trim)
  #mapview::mapview(pts_new)
  
  # get lines 
  lines <- pts_new %>%
    summarise(do_union = F) %>%
    st_cast("LINESTRING")
  
  # make plot to save and check after
  p <- ggplot()+
    geom_sf(data = iso20) +
    geom_sf(data = lines, color = "gray29", size = 0.5) +
    geom_sf(data = pts_new, color = "red", size = 0.7)
  p
  
  # ggsave(plot = p, path = here("outputs","pathroutr_plots"),
  #        filename = paste0(id,"_behavlog_pseudotrack_rerouted_20mIso_2025Feb18.png"))
  
  # remove the "original" timestamps of the argos locations, so we only get the 
  # 30min records
  id_pseu4 <- filter(pts_new, locType == "p")
  id_pseu4$DeployID <- id
  
  # store
  pseus[[id]] <- id_pseu4

  
}

# get them each individually 
pair1 <- pseus[["PcTag030"]]
pair2 <- pseus[["PcTag032"]]

## calculate pairwise distances ## ------------------------------------------- ##
pair1$dist_m <- NA

# first record of PcTag030/090 was before that of the other individual in the pair
# so remove that record
pair1 <- pair1[c(2:nrow(pair1)),]

for(i in 1:nrow(pair1)){
  #i = 1
  
  pair1[i,"dist_m"] <- as.numeric(st_distance(pair1[i,],pair2[i,], by_element = F))
  
}

# clean up the data and save for visualizations 
colnames(pair1)
pair_final <- pair1 %>%
  select(DeployID,timestamp,
         se.mu.x, se.mu.y, dist_m, geometry) %>%
  mutate(
    dist_km = dist_m/1000
  )

# save pair file and both tracks 
save(pair, pair_final, pair1, pair2, file = here("pipeline","crawl_pseudotracks",
                                           "PcTag030_PcTag032_pairwise_dists_locs_2025Apr29.RData"))

save(pair, pair_final, pair1, pair2, file = here("pipeline","crawl_pseudotracks",
                                                 "PcTag090_PcTag092_pairwise_dists_locs_2025Apr29.RData"))

# get simple summary stats to report 
median(pair_final$dist_km)
mean(pair_final$dist_km)
max(pair_final$dist_km)

# make a subset with locations with more > 1km error removed 
pair_final$avgse <- (pair_final$se.mu.x + pair_final$se.mu.y)/2
pair_sub <- filter(pair_final, avgse < 1000)

# recalculate summary stats
median(pair_sub$dist_km)
mean(pair_sub$dist_km)
max(pair_sub$dist_km)