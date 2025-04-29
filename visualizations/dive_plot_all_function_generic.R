DivePlotAll = function(behavior,
                       title = NULL,
                       gapcolor = "black",
                       plotcolor = "black",
                       filetz = "UTC",
                       plottz = "Pacific/Honolulu",
                       depths = c(-750 -500, -250, 0),
                       depthlab = c("1500", "1000", "-500", "-250", "0"),
                       depthlim = NULL,
                       datelim = NULL,
                       nights = NULL,
                       datebreaks = "2 days",
                       lineweight = 0.2,
                       gaps = T,
                       grid = T) {
  
  
  
  behavior$Start = with_tz(as.POSIXct(gsub("\\.5", "", behavior$Start), tz=filetz, format="%H:%M:%S %d-%b-%Y"), tzone=plottz)
  behavior$End = with_tz(as.POSIXct(gsub("\\.5", "", behavior$End), tz=filetz, format="%H:%M:%S %d-%b-%Y"), tzone=plottz)
  behavior$DepthMin = 0 - behavior$DepthMin
  messages = subset(behavior, What == "Message")
  surface = subset(behavior, What == "Surface")
  dive = subset(behavior, What == "Dive")
  
  p = ggplot(data = behavior, aes(Start, DepthMin)) + theme_bw() +
    theme(#axis.title.x = element_blank(), 
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
    coord_cartesian(xlim = datelim)
  if (!grid) {
    p = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }
  
  if (!is.null(nights)) {
    if (gaps) {
      p = p + geom_rect(data = nights, mapping = aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = 0),
                        fill = "grey", alpha = 0.5, inherit.aes = F)
    } else {
      p = p + geom_rect(data = nights, mapping = aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                        fill = "grey", alpha = 0.5, inherit.aes = F)
    }
  }
  
  if (!is.null(title)) {
    p = p + ggtitle(title)
  }
  
  datetime = with_tz(as.POSIXct(NA), tzone=plottz)
  depth = as.numeric(NA)
  grp = as.numeric(NA)
  
  ### NEED TO BE MORE EXPLICIT WITH TIME ZONES AND CONCATENATE.
  ## If you do not tell c() what time zone each element should be in,
  ## it will default to that on on your machine. In most cases this isn't an
  ## issue because we've already set the time zone for the Start column
  ## However as soon as you convert a duration to a numeric and attempt to add
  ## things get panicky. Be explicit in each case and all is Hawaii time
  for (i in 1:nrow(dive)) {
    DiveStart = dive[i, "Start"]
    DiveDur = as.numeric(as.duration(int_diff(c(dive[i, "Start"], dive[i, "End"]))))
    Depth = dive[i, "DepthMin"]
    if (dive[i, "Shape"] == "V") {
      ## ADDED WITH_TZ
      datetime = with_tz(c(datetime, DiveStart, DiveStart + (DiveDur * 0.5), DiveStart + DiveDur), tzone=plottz)
      depth = c(depth, 0, Depth, 0)
      grp = c(grp, i, i, i)
    } else if (dive[i, "Shape"] == "U") {
      ## ADDED WITH_TZ
      datetime = with_tz(c(datetime, DiveStart, DiveStart + (DiveDur * 0.25), DiveStart + (DiveDur * 0.75), DiveStart + DiveDur), tzone=plottz)
      depth = c(depth, 0, Depth, Depth, 0)
      grp = c(grp, i, i, i, i)
    } else if (dive[i, "Shape"] == "Square") {
      ## ADDED WITH_TZ
      datetime = with_tz(c(datetime, DiveStart, DiveStart + (DiveDur * 0.4), DiveStart + (DiveDur * 0.6), DiveStart + DiveDur), tzone=plottz)
      depth = c(depth, 0, Depth, Depth, 0)
      grp = c(grp, i, i, i, i)
    }
  }
  
  dive.df = data.frame(datetime, depth, grp)
  p = p + geom_line(data = dive.df, mapping = aes(datetime, depth, group =  grp), linewidth = lineweight, color = plotcolor)
  
  datetime = with_tz(as.POSIXct(NA), tzone=plottz)
  depth = as.numeric(NA)
  grp = as.numeric(NA)
  
  for (i in 1:nrow(surface)) {
    datetime = with_tz(c(datetime, surface[i, "Start"], surface[i, "End"]), tzone=plottz)
    depth = c(depth, 0, 0)
    grp = c(grp, i, i)
  }
  
  surface.df = data.frame(datetime, depth, grp)
  p = p + geom_line(data = surface.df, mapping = aes(datetime, depth, group =  grp), linewidth = lineweight, color = plotcolor)
  
  if (is.null(depthlim)) {
    depth_lower = min(dive$DepthMin) * 1.15
  } else {
    depth_lower = depthlim
  }
  if (gaps) {
    depth_upper = abs(depth_lower) * 0.15
  } else {
    depth_upper = 0
  }
  
  xmin = with_tz(as.POSIXct(NA), tzone=plottz)
  xmax = with_tz(as.POSIXct(NA), tzone=plottz)
  
  gapCheck = 0
  for (i in 2:nrow(messages)) {
    
    dt <- difftime( messages[i-1, "End"], messages[i, "Start"], units = 'secs')
    ## Need an indicator that avoids plotting gaps if none exist
    ## Otherwise you get screwed with gaps = TRUE
    if (!is.na(dt)) {
      if (dt < -60) {
        xmin = c(xmin, messages[i-1, "End"])
        xmax = c(xmax, messages[i, "Start"])
        ymin = abs(depth_lower) * 0.05
        ymax = abs(depth_lower) * 0.15
        ## Update if gaps were found
        gapCheck = gapCheck + 1
      }
    }
  }
  
  gapCheck <- ifelse(gapCheck > 0, 1, 0) 
  
  
  ## If you are plotting gaps
  if (gaps & gapCheck) {
    gaps.df = data.frame(xmin, xmax, ymin, ymax)
    gaps.df = gaps.df[!is.na(gaps.df$xmin), ]
    p = p + geom_rect(data = gaps.df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                      fill = gapcolor, inherit.aes = F) +
      scale_y_continuous(breaks = c(depths, mean(c(ymin, ymax))),
                         labels = c(depthlab, "Gaps"),
                         limits = c(depth_lower, depth_upper)) +
      scale_x_datetime(labels = date_format("%Y-%m-%d", tz = plottz), breaks = date_breaks(datebreaks),
                       expand = c(0,0)) +
      ylab("Depth (m)")
  } else {
    p = p + scale_y_continuous(breaks = depths,
                               labels = depthlab,
                               limits = c(depth_lower, depth_upper)) +
      scale_x_datetime(labels = date_format("%Y-%m-%d", tz = plottz), breaks = date_breaks(datebreaks),
                       expand = c(0,0)) +
      ylab("Depth (m)")
  }
  return(p)
}
