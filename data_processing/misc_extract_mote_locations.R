## extract_mote_locations.R: extract mote locations that SPLASH-tagged FKWs
## had hits from for mapping 

## Author: Michaela A. Kratofil, Oregon State University, Cascadia Research
## Updated: 21 Feb 2025

## --------------------------------------------------------------------------- ##

## load packages 
library(dplyr)
library(lubridate)
library(here)

## read in -all.csv files for tags that were within range of islands ## ------ ##
files <- list.files(path = here("data","location"),pattern = "-All.csv")
files

# for loop to read in files 
motes_list <- list()
for(i in 1:length(files)){
  # for testing 
  #i = 1
  
  # read in the file
  f <- read.csv(here("data","location", files[i]))
  
  # get the deploy id
  dep <- f$DeployID[1]
  
  # get list of unique mote ids 
  mote_messages <- filter(f, Loc..type == "Mote")
  
  # create dataframe for mote messages, even if there were none
  if(nrow(mote_messages) > 0){
 
    mote_df <- mote_messages %>%
      group_by(Mote.Id) %>%
      slice(1) %>%
      select(DeployID, Mote.Id, Latitude, Longitude)
  } else{
    
    mote_df <- data.frame(DeployID = dep,
                          Mote.Id = NA,
                          Latitude = NA,
                          Longitude = NA)
    
  }
  
  # store the dataframe in the list 
  motes_list[[i]] <- mote_df

}

# bind list 
motes <- bind_rows(motes_list)

# save the list 
write.csv(motes, here("pipeline","fkw_splash_mote_ids_by_deploy_id.csv"), row.names = F)
