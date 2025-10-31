## pseudotrack.R: Fit SPLASH Pseudorca tag data to crawl and predict locations
## at each time stamp of the behavior log, effectively creating a "pseudotrack"

## 02a: regular/single pseudotrack (not multiple imputation)

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 04 Oct 2025

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

# set seed 
set.seed(123)

## read in and prep Douglas Filtered movement data for modeling ## ----------- ##
tbl_locs <- readr::read_csv(here("data","location","PcTag001-099_DouglasFiltered_r20d3lc2_ArgosGPS_2025OCTv1.csv"),
                            col_types = 
                            cols(animal = col_character(),
                            ptt = col_double(),
                            date = col_datetime(),
                            longitud = col_double(),
                            latitude = col_double(),
                            LC = col_character(),
                            error_radius = col_double(),
                            ellipse_orient = col_double(),
                            semi_major = col_double(),
                            semi_minor = col_double(),
                            sat_count = col_integer(),
                            LocType = col_character()))
tbl_locs

# review
summary(tbl_locs)
length(unique(tbl_locs$animal))
class(tbl_locs$date)
tz(tbl_locs$date)

# check lcs
unique(tbl_locs$LC)

# select SPLASH tags only 
tbl_sub <- tbl_locs %>%
  filter(animal %in% c("PcTag026","PcTag028","PcTag030","PcTag032",
                       "PcTag035","PcTag037","PcTag049","PcTag055",
                       "PcTag074","PcTag090","PcTag092","PcTagP09",
                       "PcTag095","PcTag096","PcTag097","PcTag099")) 

summary(tbl_sub)

# clean up data, convert everything to lower case and rename columns
tbl_sub <- tbl_sub %>%
  janitor::clean_names() %>% # this will drop everything to lower case to avoid typos
  dplyr::rename(lon = longitud,
                lat = latitude,
                timestamp = date,
                e_radius = error_radius,
                e_orient = ellipse_orient,
                smaj = semi_major,
                smin = semi_minor) %>%
  dplyr::arrange(animal, timestamp)
head(tbl_sub)

# set error radius to 50 m for deployment (DP) locations, known from GPS coordinates from boat
unique(tbl_sub$sat_count)
tbl_sub <- tbl_sub %>%
  mutate(
    smaj = ifelse(lc == "DP", 50, smaj),
    smin = ifelse(lc == "DP", 50, smin),
    e_orient = ifelse(lc %in% c("DP","G"), 0, e_orient),
    smaj = ifelse(lc == "G" & sat_count > 4, 50, smaj),
    smaj = ifelse(lc == "G" & sat_count == 4, 1000, smaj),
    smin = ifelse(lc == "G" & sat_count > 4, 1000, smin),
    smin = ifelse(lc == "G" & sat_count == 4, 1000, smin)
  )

# check
summary(tbl_sub)

# remove location without error ellipse info
tbl_sub <- filter(tbl_sub, !is.na(e_orient))

# Project data prior to fitting to crawl model. 
sf_locs <- sf::st_as_sf(tbl_sub, coords = c("lon","lat")) %>% # spatial points 
  sf::st_set_crs(4326) %>%
  sf::st_transform(crs = 3857)

# set fixed params for all models. First 2 will be 1 (location uncertainty), and sigma and beta will be NA
# since we'll be estimating those 
fixPar <- c(1,1,NA,NA)

# nest data by tag for tidy fitting
mov <- sf_locs %>% 
  dplyr::group_by(animal, ptt) %>% 
  tidyr::nest() %>% ungroup()

# create functions to identify periods within deployments with large gaps 
# between argos data 
get_segments <- function(x, gap = 4, time_unit = "hours"){
  dt <- diff(x) %>% `units<-`(time_unit)
  time_diff <- c(0, dt)
  seg_id <- (time_diff >= gap) %>% {cumsum(.)+1}
  return(seg_id)
}

location_rate <- function(x, time_unit="hour", stat=mean, ...){
  out <- table(lubridate::round_date(x$timestamp, time_unit))
  return(stat(out,...))
}

# identify 4, 6, and 8 hr gaps in the data
mov <- mov %>% mutate(
  data = map(data, ~{
    .x$seg_id_4 = get_segments(.x$timestamp,4)
    .x$seg_id_6 = get_segments(.x$timestamp,6)
    .x$seg_id_8 = get_segments(.x$timestamp,8)
    .x }),
  num_seg_4 = map_dbl(data, ~{max(.x$seg_id_4)}),
  num_seg_6 = map_dbl(data, ~{max(.x$seg_id_6)}),
  num_seg_8 = map_dbl(data, ~{max(.x$seg_id_8)})
)

# create location uncertainty matrix to feed into crawl
mov <- mov %>% 
  dplyr::mutate(
    diag = purrr::map(data, ~ crawl::argosDiag2Cov(
      .x$smaj, 
      .x$smin, 
      .x$e_orient)),
    data = purrr::map2(data,diag,bind_cols)) %>% 
  dplyr::select(-diag)


# fit CTCRW movement model (same set up as what Devin runs on FKW data) ## --- ##
fit1_func <- function(d, fixPar) {
  suppressWarnings(
    crwMLE(mov.model = ~ 1,
           err.model = list(
             x = ~ ln.sd.x + 0,
             y = ~ ln.sd.y + 0,
             rho = ~ error.corr
           ),
           fixPar = fixPar,
           data = d,
           theta = c(9,0),
           Time.name = 'timestamp',
           attempts = 10,
           control = list(maxit = 2000, trace = 0, REPORT = 10),
           initialSANN = list(maxit = 1500, trace = 0, REPORT = 10, temp = 10))
  )
}

# apply the function (model for each individual)
mov <- mov %>%
  dplyr::mutate(fit1 = purrr::map(data, ~ fit1_func(d = .x, fixPar = fixPar)))

# check individual parameter estimates (can print them out by switching the
# second indexed value from 1 (first tag), 2 (second tag), etc.)
for(i in 1:nrow(mov)){
  print(mov[[7]][[i]])
}

# all params look good, use this model.
# save crawl fitted models RDS file
saveRDS(mov, here("pipeline", "PcTags_SPLASH_thru099_Kalman_Crawl_Fitted_Model_wSegments_2025Oct04.rds"))

## predict crawl positions at behavior log times ## -------------------------- ##
mov <- readRDS(here("pipeline","PcTags_SPLASH_thru099_Kalman_Crawl_Fitted_Model_wSegments_2025Oct04.rds"))
mov

# get objects for pathroutr to reroute predicted locations around land 
# read in 20m isobath 
iso20 <- readRDS(here("code","data_processing","isobath20m_for_rerouting.rds"))
bathy <- readRDS(here("code","data_processing","raster_for_rerouting.rds"))

# quick map
ggplot()+
  geom_sf(data=iso20)+
  theme_classic()+ theme(legend.position="none")
st_crs(iso20) # check crs

# pathroutr: visability graph for rerouting paths around barrier (iso20)
iso20_b <- st_buffer(iso20, 10000)

# pathroutr: add augmented points for isobaths that aren't as complete, like niihau,
#  kauai, nihoa, and necker
aug_pts <- raster::rasterToContour(bathy, levels=c(-100, -200, -300)) %>% st_as_sf() %>% 
  st_sample(300, type='regular') %>% st_cast("POINT")
aug_pts <- aug_pts[st_intersects(aug_pts, iso20_b, sparse=F)]  

# quick check of aug points 
ggplot() +
  geom_sf(data = iso20) +
  geom_sf(data = aug_pts)

# pathroutr: create visibility graph
vis_graph <- pathroutr::prt_visgraph(iso20, centroids=TRUE, aug_points=aug_pts)

# going to loop through each animal, read in the behavior log file, and provide
# the timestamps of the behavior log file to crawl to predict locations at. will
# also predict at finer timestep (this helps ensure the pathroutr process works
# as intended)
ids <- unique(mov$animal)

for(i in 1:length(ids)){
  # for testing
  #i = 3
  
  # select the animal
  id <- ids[i]
  
  # get the compiled QA/QC'd behavior log file 
  bq <- readRDS(here("pipeline","qaqcd_behavior_data",
                     "PcTags_SPLASH_thru099_behavior_logs_qaqcd_2025Oct02.rds"))
  
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
  
  # create a vector of times every 5 mins from the start to the end of the 
  # behavior log
  times_pred <- data.frame(start_utc = timeplyr::time_seq(id_dive$start_utc[1],
                                                         id_dive$start_utc[nrow(id_dive)],
                                                         time_by = "5 min"),
                           type = "pred") %>%
    bind_rows(., id_dive) %>%
    arrange(start_utc) %>%
    mutate(
      dt = as.numeric(difftime(start_utc, lag(start_utc), units = "mins"))
    ) %>%
    dplyr::select(start_utc, type, dt) 
  
  # filter duplicates 
  nodup <- times_pred %>%
    filter(., dt > 0 | is.na(dt))
  
  # create vector of times to predict crawl for
  log_times <- nodup$start_utc
  
  # get the crawl fitted object for that ID
  id_crawl <- mov[[7]][[i]]
  
  # predict crawl locations for the behavior log times 
  id_pseu <- crwPredict(id_crawl, predTime = log_times, return.type = "flat") %>%
    crw_as_sf(.,  ftype = "POINT", locType = c("p","o"))
  id_pseu$start_utc <- id_pseu$timestamp
  #mapview::mapview(id_pseu)
  
  # append the other variables from dive behavior log to the crawl file 
  id_pseu2 <- left_join(id_pseu, id_dive, by = "start_utc")
  
  # identify segments (gaps between big argos locations) now with the behavior
  # log data appended
  id_pseu3 <- id_pseu2 %>%
    mutate(
      tmp_seg4.1 = zoo::na.locf(seg_id_4),
      tmp_seg4.2 = zoo::na.locf(seg_id_4, fromLast = T),
      seg_id4 = ifelse(tmp_seg4.1 != tmp_seg4.2, NA, tmp_seg4.1),
      
      tmp_seg6.1 = zoo::na.locf(seg_id_6),
      tmp_seg6.2 = zoo::na.locf(seg_id_6, fromLast = T),
      seg_id6 = ifelse(tmp_seg6.1 != tmp_seg6.2, NA, tmp_seg6.1),
      
      tmp_seg8.1 = zoo::na.locf(seg_id_8),
      tmp_seg8.2 = zoo::na.locf(seg_id_8, fromLast = T),
      seg_id8 = ifelse(tmp_seg8.1 != tmp_seg8.2, NA, tmp_seg8.1)
    ) 
  
  # reroute the locations around land 
  pts_trim <- pathroutr::prt_trim(id_pseu3, iso20)
  pts_rr <- pathroutr::prt_reroute(pts_trim, iso20, vis_graph, blend = F)
  
  # manual workaround bug in pathrotur for prt_update:
  # if there were rerouted points, then update geometry. otherwise, use pts trim
  if(nrow(pts_rr) > 0){
    
    # update geometry of points with rerouted geometry. doing this manually for now
    # because keep getting an error with dplyr::rows_update() (told Josh)
    
    # first get a fixed/point id of the trimmed points
    pts_trim <- pts_trim %>% sf::st_cast("POINT") %>%
      sf::st_sf() %>% tibble::rowid_to_column("fid")
    
    # get st coordinates and convert to dataframe
    trim_coords <- st_coordinates(pts_trim)
    trim_df <- pts_trim %>%
      st_drop_geometry() %>%
      bind_cols(., trim_coords)
    trim_df$or <- "original" # indicator for type (helpful for checking later)
    
    # rerouted points: not in sf object, so make into sf object, then dataframe
    rr_sf <- pts_rr %>%
      st_sf(sf_column_name = "geometry", crs = st_crs(pts_trim))
    rr_coords <- st_coordinates(rr_sf)
    rr_df <- rr_sf %>%
      st_drop_geometry() %>%
      bind_cols(., rr_coords)
    rr_df$or <- "rerouted" # indicator for type (helpful for checking later)
    
    # update the rows in trimmed dataframe with rerouted points, based on fid
    res <- trim_df %>%
      dplyr::rows_update(rr_df, by = "fid")
    
    # make into sf object again
    res_sf <- st_as_sf(res, coords = c("X","Y"), crs = st_crs(pts_trim))
    
    # check 
    # ggplot() +
    #   geom_sf(data = iso20, color = NA, fill = "black") +
    #   geom_sf(data = pts_trim, color = "red") +
    #   geom_sf(data = res_sf, color = "blue")
    
    
  } else{
    
    # first get a fixed/point id of the trimmed points
    pts_trim <- pts_trim %>% sf::st_cast("POINT") %>%
      sf::st_sf() %>% tibble::rowid_to_column("fid")
    
    # get st coordinates and convert to dataframe
    trim_coords <- st_coordinates(pts_trim)
    trim_df <- pts_trim %>%
      st_drop_geometry() %>%
      bind_cols(., trim_coords)
    trim_df$or <- "original" # indicator for type (helpful for checking later)
    
    # make into sf object again (same format as above)
    res_sf <- st_as_sf(trim_df, coords = c("X","Y"), crs = st_crs(pts_trim))
    
  }
  
  # reassign
  pts_new <- res_sf
  
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
  
  ggsave(plot = p, path = here("outputs","pathroutr_plots"),
         filename = paste0(id,"_behavlog_pseudotrack_rerouted_20mIso_2025Oct04.png"))
  
  # remove the "original" timestamps of the argos locations, so we only get the 
  # dive/surface records
  id_pseu4 <- filter(pts_new, !is.na(What))
  
  # get the coordinates back
  id_pseu4_wgs <- st_transform(id_pseu4, crs = 4326)
  id_coords <- as.data.frame(st_coordinates(id_pseu4_wgs))
  
  # format/clean up the data 
  clean_pseu <- id_pseu4 %>%
    bind_cols(., id_coords) %>%
    st_drop_geometry() %>%
    dplyr::rename(
      lon = X,
      lat = Y
    ) %>%
    dplyr::select(DeployID, Ptt, timestamp, Start, End, What, Shape, DepthMin, DepthMax,
                  DurationMin, DurationMax, Shallow, Deep, lon, lat, se.mu.x, se.mu.y, seg_id4,
                  seg_id6, seg_id8, start_utc, end_utc, start_hst) %>%
    rename(
      crawl_timestamp = timestamp
    )
  
  # save the file 
  write.csv(clean_pseu, here("pipeline","crawl_pseudotracks",
                             paste0(id,"_behavlog_pseudotrack_rerouted20mIso_2025Oct04.csv")), row.names = F)
  
}
