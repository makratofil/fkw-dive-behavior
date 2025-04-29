split_surfs_tod <- function(surfs, bh){
  
  for(i in 1:nrow(surfs)) {
    cat("\r", i, "                     ")
    
    # for testing 
    #surfs <- surfs_nest[[4]][[1]]
    #bh <- surfs_nest[[5]][[1]]
    #i = 1
    
    # get basics info
    DeployID = surfs[i, "DeployID", drop = T]
    Ptt = surfs[i, "Ptt", drop = T]
    TZ = surfs[i, "TZ", drop = T]
    
    # format start and end times for this tag 
    #surfs$Start <- as.POSIXct(surfs$Start, tz = TZ)
    #surfs$End <- as.POSIXct(surfs$End, tz = TZ)
    
    # bind empty dataframe with ith row of surface data
    bh = bind_rows(bh, surfs[i, ])
    
    # get corresponding times for ith row of surface data
    dt0 = surfs[i, "datetime_hst", drop = T]
    dt1 = as.POSIXct(format(surfs[i, "End", drop = T], tz="Pacific/Honolulu"), tz="Pacific/Honolulu")
    da_st = surfs[i, "civil_dawn", drop = T]
    da_en = surfs[i, "end_dawn", drop = T]
    du_st = surfs[i, "start_dusk", drop = T]
    du_en = surfs[i, "civil_dusk", drop = T]
    
    # get the number of days that the surface interval lasted 
    surf_dur <- ceiling(surfs[i, "durationDays", drop = T])
    
    # create dataframe for the dawn and dusk times for the day of the surface
    # period ("day0") and subsequent days of the surface interval
    surf_days <- data.frame(day = c(0, rep(NA, surf_dur)),
                            civil_dawn = c(da_st, rep(NA, surf_dur)),
                            end_dawn = c(da_en, rep(NA, surf_dur)),
                            start_dusk = c(du_st, rep(NA, surf_dur)),
                            civil_dusk = c(du_en, rep(NA, surf_dur)))
    
    # get dawn and dusk times for the sequence of days that the surface period
    # covers IF it is more than one day long 
    #if(surf_dur > 1) {
    for(j in 2:nrow(surf_days)){
      #j = 2
      surf_days[j, "day"] <- surf_days[j-1, "day"] + 1
      surf_days[j, "civil_dawn"] <- surf_days[j-1, "civil_dawn"] + days(1)
      surf_days[j, "end_dawn"] <- surf_days[j-1, "end_dawn"] + days(1)
      surf_days[j, "start_dusk"] <- surf_days[j-1, "start_dusk"] + days(1)
      surf_days[j, "civil_dusk"] <- surf_days[j-1, "civil_dusk"] + days(1)
      
    }
    #}
    
    # create intervals for the surf days
    surf_days <- surf_days %>%
      mutate(
        dawn_int = interval(civil_dawn, end_dawn, tzone = tz(dt0)),
        dusk_int = interval(start_dusk, civil_dusk, tzone = tz(dt0))
        
      )
    
    # pivot dataframe longer so that it's easier to identify which of the
    # dawn/dusk intervals the surface interval overlaps with 
    surf_days_long <- surf_days %>%
      tidyr::pivot_longer(cols = c(dawn_int, dusk_int), names_to = "type",
                          values_to = "int") %>%
      mutate(
        dday = paste0(type, day),
        tod = case_when(stringr::str_detect(dday, "dawn") ~ "dawn", TRUE ~ "dusk"),
        start = case_when(stringr::str_detect(dday, "dawn") ~ civil_dawn, TRUE ~ start_dusk),
        end = case_when(stringr::str_detect(dday, "dawn") ~ end_dawn, TRUE ~ civil_dusk)
      ) 
    
    # create interval for ith surface interval
    surf_int <- interval(start = dt0, end = dt1, tzone = tz(dt0))
    
    # now determine which of the dawn/dusk intervals the surface interval overlaps
    # with 
    surf_days_long <- surf_days_long %>%
      mutate(
        overlap = int_overlaps(surf_int, int)
      )
    
    # if the surface interval overlaps with any of the dawn/dusk intervals, then
    # move through the splitting process. if it does not overlap with any
    # dawn/dusk intervals, then assign time of day and move to the next surface 
    # interval 
    if(any(surf_days_long$overlap == T, na.rm = T)){
      
      # get the intervals that that surface interval overlaps with
      overlap_ints <- filter(surf_days_long, overlap == T)
      
      # determine whether the start of the surface interval occurs before or 
      # during the first dawn/dusk interval, and whether the end of the surface 
      # interval ends during or after the last dawn/dusk interval 
      
      # OPTION 1: starts before the first dawn/dusk interval & ends before
      # the end of the n-th dawn/dusk interval 
      if(dt0 < first(overlap_ints$start) && dt1 < last(overlap_ints$end)){
        
        # assign the vector of datetimes to make the splits: in this case, this
        # will be the start of the surface interval, all start and end of dawn/
        # dusk times after EXCEPT for the last end time, and the end of the
        # surface interval. if there is only one overlapping interval, then this
        # will simply be dt1
        if(nrow(overlap_ints) == 1){
          int_dates_df <- data.frame(dates = c(dt0, overlap_ints$start, dt1)) %>%
            arrange(dates) # need to arrange them in order of time first
          
          int_dates <- int_dates_df$dates # now create vector to provide in function
        } else{
          int_dates_df <- data.frame(dates = c(dt0, overlap_ints$start, overlap_ints$end[1:nrow(overlap_ints)-1], dt1)) %>%
            arrange(dates) # need to arrange them in order of time first
          
          int_dates <- int_dates_df$dates # now create vector to provide in function
        }
        # determine the sequence of "time of day" to assign to each split. first
        # determine the start and end time of day. 
        if(first(overlap_ints$tod) == "dawn"){
          dt0_tod <- "night"
        } else{
          dt0_tod <- "day"
        }
        
        if(last(overlap_ints$tod) == "dawn"){
          dt1_tod <- "dawn"
        } else{
          dt1_tod <- "dusk"
        }
        
        # now we have to get the sequence of "dawn, day, dusk, night" in between
        # the tod of the start and end of the surface interval. 
        # add NA rows to every other row in overlap_ints so we can add "day" and
        # "night" in between the dawn and dusk periods 
        tods_df <- overlap_ints[rep(1:nrow(overlap_ints), each = 2), ]
        tods_df[1:nrow(tods_df) %% 2 == 0, ] <- NA
        
        # because this interval ends before the end of the dawn/dusk period,
        # we'll want to remove the last added row
        tods_df <- tods_df[1:nrow(tods_df)-1,]
        
        # now fill in the missing rows with day and night
        tods_df <- tods_df %>%
          mutate(
            tod = ifelse(is.na(tod) & lag(tod) == "dusk", "night", tod),
            tod = ifelse(is.na(tod) & lag(tod) == "dawn", "day", tod)
          )
        
        # now create the entire vector of all tods. in this case, the operation
        # above created the last record (dt1_tod), so we don't need to append
        # dt1_tod
        tods_split <- c(dt0_tod, tods_df$tod)
        
        # split up periods with split surfaces function
        bh <- split_surfs(bh, dates = int_dates, tods = tods_split, tag_tz = TZ,
                          dep_id = DeployID, ptt = Ptt)
        
        
        
      } 
      # OPTION 2: starts before the first dawn/dusk interval & ends after the 
      # end of the n-th dawn/dusk interval
      else if(dt0 < first(overlap_ints$start) && dt1 > last(overlap_ints$end)){
        
        # assign the vector of datetimes to make the splits: in this case, this
        # will be the start of the surface interval, all start and end of dawn/
        # dusk times after, and the end of the surface interval.  
        int_dates_df <- data.frame(dates = c(dt0, overlap_ints$start, overlap_ints$end, dt1)) %>%
          arrange(dates) # need to arrange them in order of time first
        
        int_dates <- int_dates_df$dates # now create vector to provide in function
        
        # determine the sequence of "time of day" to assign to each split. first
        # determine the start and end time of day. 
        if(first(overlap_ints$tod) == "dawn"){
          dt0_tod <- "night"
        } else{
          dt0_tod <- "day"
        }
        
        if(last(overlap_ints$tod) == "dawn"){
          dt1_tod <- "day"
        } else{
          dt1_tod <- "night"
        }
        
        # now we have to get the sequence of "dawn, day, dusk, night" in between
        # the tod of the start and end of the surface interval. 
        # add NA rows to every other row in overlap_ints so we can add "day" and
        # "night" in between the dawn and dusk periods 
        tods_df <- overlap_ints[rep(1:nrow(overlap_ints), each = 2), ]
        tods_df[1:nrow(tods_df) %% 2 == 0, ] <- NA
        
        # now fill in the rows with day and night
        tods_df <- tods_df %>%
          mutate(
            tod = ifelse(is.na(tod) & lag(tod) == "dusk", "night", tod),
            tod = ifelse(is.na(tod) & lag(tod) == "dawn", "day", tod)
          )
        
        # now create the entire vector of all tods. in this case, the operation
        # above created the last record (dt1_tod), so we don't need to append
        # dt1_tod
        tods_split <- c(dt0_tod, tods_df$tod)
        
        # split up periods with split surfaces function
        bh <- split_surfs(bh, dates = int_dates, tods = tods_split, tag_tz = TZ,
                          dep_id = DeployID, ptt = Ptt)
        
      }
      # OPTION 3: starts after the first dawn/dusk interval & ends before the 
      # end of the n-th dawn/dusk interval
      else if(dt0 > first(overlap_ints$start) && dt1 < last(overlap_ints$end)){
        
        # assign the vector of datetimes to make the splits: in this case, this
        # will be the start of the surface interval, all start and end of dawn/
        # dusk times after EXCEPT for the first dawn/dusk start time AND the
        # last dawn/dusk end time, and the end of the surface interval. if there
        # is only one overlapping interval, then this will be dt0 and dt1 and
        # nothing needs to be split.
        if(nrow(overlap_ints) == 1){
          
          # no splitting, just assign time of day 
          bhrow = nrow(bh)
          bh[bhrow, "tod"] <- overlap_ints$tod
          
          # add columns for start and end in HST just like we did for the splits 
          bh[bhrow, "start_hst"] <- as.POSIXct(format(bh[bhrow, "Start", drop = T], tz = "Pacific/Honolulu"), tz = "Pacific/Honolulu")
          bh[bhrow, "end_hst"] <- as.POSIXct(format(bh[bhrow, "End", drop = T], tz = "Pacific/Honolulu"), tz = "Pacific/Honolulu")
          
        } else{
          int_dates_df <- data.frame(dates = c(dt0, overlap_ints$start[2:nrow(overlap_ints)],
                                               overlap_ints$end[1:nrow(overlap_ints)-1], dt1)) %>%
            arrange(dates) # need to arrange them in order of time first
          
          int_dates <- int_dates_df$dates # now create vector to provide in function
          
          # determine the sequence of "time of day" to assign to each split. first
          # determine the start and end time of day. IN THIS CASE: the start TOD
          # is equal to the first dawn/dusk interval, and the last TOD is equal
          # to the last dawn/dusk interval (surface period is completely within
          # the respective dawn/dusk intervals)
          
          # add NA rows to every other row in overlap_ints so we can add "day" and
          # "night" in between the dawn and dusk periods 
          tods_df <- overlap_ints[rep(1:nrow(overlap_ints), each = 2), ]
          tods_df[1:nrow(tods_df) %% 2 == 0, ] <- NA
          
          # remove the last row, because the surface interval ends before the next
          # TOD in the sequence
          tods_df <- tods_df[1:nrow(tods_df)-1,]
          
          # now fill in the NA rows with day and night 
          tods_df <- tods_df %>%
            mutate(
              tod = ifelse(is.na(tod) & lag(tod) == "dusk", "night", tod),
              tod = ifelse(is.na(tod) & lag(tod) == "dawn", "day", tod)
            )
          
          # now create the entire vector of all tods. in this case, it will only
          # be the existing TODs
          tods_split <- c(tods_df$tod)
          
          # split up periods with split surfaces function
          bh <- split_surfs(bh, dates = int_dates, tods = tods_split, tag_tz = TZ,
                            dep_id = DeployID, ptt = Ptt)
        }
        
        
        
      }
      # OPTION 4: starts after the first dawn/dusk interval & ends after the 
      # end of the n-th dawn/dusk interval
      else if(dt0 > first(overlap_ints$start) && dt1 > last(overlap_ints$end)){
        
        # assign the vector of datetimes to make the splits: in this case, this
        # will be the start of the surface interval, all start and end of dawn/
        # dusk times after EXCEPT for the first dawn/dusk start time, and the 
        # end of the surface interval. if there is only one overlapping interval,
        # then this will be the start time, end of the dawn/dusk period, and end
        # time
        if(nrow(overlap_ints) == 1){
          
          int_dates_df <- data.frame(dates = c(dt0, overlap_ints$end, dt1)) %>%
            arrange(dates) # need to arrange them in order of time first
          
          int_dates <- int_dates_df$dates # now create vector to provide in function
          
        } else{
          int_dates_df <- data.frame(dates = c(dt0, overlap_ints$start[2:nrow(overlap_ints)], overlap_ints$end, dt1)) %>%
            arrange(dates) # need to arrange them in order of time first
          
          int_dates <- int_dates_df$dates # now create vector to provide in function
        }
        
        # determine the sequence of "time of day" to assign to each split. first
        # determine the start and end time of day. 
        if(first(overlap_ints$tod) == "dawn"){
          dt0_tod <- "dawn"
        } else{
          dt0_tod <- "dusk"
        }
        
        if(last(overlap_ints$tod) == "dawn"){
          dt1_tod <- "day"
        } else{
          dt1_tod <- "night"
        }
        
        # add NA rows to every other row in overlap_ints so we can add "day" and
        # "night" in between the dawn and dusk periods 
        tods_df <- overlap_ints[rep(1:nrow(overlap_ints), each = 2), ]
        tods_df[1:nrow(tods_df) %% 2 == 0, ] <- NA
        
        # now fill in the NA rows with day and night 
        tods_df <- tods_df %>%
          mutate(
            tod = ifelse(is.na(tod) & lag(tod) == "dusk", "night", tod),
            tod = ifelse(is.na(tod) & lag(tod) == "dawn", "day", tod)
          )
        
        # now create the entire vector of all tods. in this case, it will only
        # be the existing TODs, because the start of the interval occurs during
        # a dawn/dusk period, and the end occurs during the TOD in the following
        # level in the sequence which was determined by the operation above 
        tods_split <- c(tods_df$tod)
        
        # split up periods with split surfaces function
        bh <- split_surfs(bh, dates = int_dates, tods = tods_split, tag_tz = TZ,
                          dep_id = DeployID, ptt = Ptt)
        
      } 
    } else {
      # surface interval doesn't overlap with dawn or dusk of start day or onwards.
      # so, assign time of day variable based on whether the surface interval overlaps
      # with day or night.
      
      # first get the ith row 
      bhrow = nrow(bh)
      
      # day interval: 1 minute after the end of dawn to 1 minute before the start of dusk
      day_int <- interval(start = (da_en + 1), end = (du_st - 1), tzone = tz(da_en))
      
      # night interval for morning: midnight to 1 minute before the start of dawn
      night_am_int <- interval(start = as.POSIXct(paste0(date(dt0), "00:00:01"),
                                                  format = "%Y-%m-%d %H:%M:%S", tz = "Pacific/Honolulu"),
                               end = (da_st - 1), tzone = tz(da_st))
      
      # night interval for evening: 1 minute after the end of dusk to 1 minute before midnight 
      night_pm_int <- interval(start = (du_en + 1), end = as.POSIXct(paste0(date(dt0), "23:59:59"),
                                                                     format = "%Y-%m-%d %H:%M:%S",
                                                                     tz = "Pacific/Honolulu"), 
                               tzone = tz(da_st))
      
      # add columns for start and end in HST just like we did for the splits 
      bh[bhrow, "start_hst"] <- as.POSIXct(format(bh[bhrow, "Start", drop = T], tz = "Pacific/Honolulu"), tz = "Pacific/Honolulu")
      bh[bhrow, "end_hst"] <- as.POSIXct(format(bh[bhrow, "End", drop = T], tz = "Pacific/Honolulu"), tz = "Pacific/Honolulu")
      
      # test surface interval and assign tod accordingly
      if(surf_int %within% day_int) {
        bh[bhrow, "tod"] <- "day"
      } else if (surf_int %within% night_am_int){
        bh[bhrow, "tod"] <- "night"
      } else if (surf_int %within% night_pm_int){
        bh[bhrow, "tod"] <- "night"
      } else {
        bh[bhrow, "tod"] <- "night" # this is for intervals crossing midnight dateline (always will be night)
      }
      
    }
    
    
  }
  return(bh)
  
}

split_surfs <- function(
    bh,
    dates = NULL,
    tods = c("dusk","night","dawn","day","dusk"),
    tag_tz = "UTC",
    dep_id = NULL,
    ptt = NULL) {
  
  # for testing 
  # dates <- int_dates
  # tods <- tods_split
  # tag_tz <- "UTC"
  # dep_id <- "PcTag090"
  # ptt <- 252367
  #bh <- bh[1:43,]
  
  
  # split the periods with those dates using int_diff
  splits <- as.data.frame(int_diff(dates))
  splits$tod <- tods
  splits$start <- int_start(splits$`int_diff(dates)`)
  splits$end <- int_end(splits$`int_diff(dates)`)
  
  # shape up 
  splits_df <- splits %>%
    rename(
      start_hst = start,
      end_hst = end
    ) %>%
    mutate(
      Start = as.POSIXct(format(start_hst, tz = tag_tz), tz = tag_tz),
      End = as.POSIXct(format(end_hst, tz = tag_tz), tz = tag_tz),
      DeployID = dep_id,
      Ptt = ptt,
      What = "Surface"
    )
  
  # change the first row of bh to the  end times of first blah
  bh[nrow(bh), "End"] <- as.POSIXct(format(splits_df[1, "end_hst", drop = T], tz = "UTC"), tz = "UTC")
  bh[nrow(bh), "tod"] <- splits_df[1, "tod", drop = T]
  bh[nrow(bh), "start_hst"] <- splits_df[1,"start_hst", drop = T]
  bh[nrow(bh), "end_hst"] <- splits_df[1,"end_hst", drop = T]
  
  # now append the rest of the splits
  bhs <- bind_rows(bh, splits_df[2:nrow(splits_df),])
  
  # return the final product 
  return(bhs)
  
  
}