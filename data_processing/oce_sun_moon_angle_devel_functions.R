sunAngle <- function(t, longitude = 0.0, latitude = 0.0, useRefraction = FALSE) {
  if (missing(t)) {
    stop("must provide t")
  } else {
    if (is.character(t)) {
      t <- as.POSIXct(t, tz = "UTC")
    }
    if (inherits(t, "Date")) {
      t <- as.POSIXct(t)
    }
    if (!inherits(t, "POSIXt")) {
      if (is.numeric(t)) {
        tref <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC") # arbitrary
        t <- t - as.numeric(tref) + tref
      } else {
        stop("t must be POSIXt or a number corresponding to POSIXt (in UTC)")
      }
    }
  }
  t <- as.POSIXct(t) # so we can get length ... FIXME: silly, I know
  nt <- length(t)
  nlongitude <- length(longitude)
  nlatitude <- length(latitude)
  nuseRefraction <- length(useRefraction)
  if (nlongitude != nlatitude) {
    stop("lengths of longitude and latitude must match")
  }
  if (nlongitude != nuseRefraction) {
    stop("lengths of longitude and useRefraction must match")
  }
  if (nlongitude == 1) {
    longitude <- rep(longitude, length.out = nt)
    latitude <- rep(latitude, length.out = nt)
    useRefraction <- rep(useRefraction, length.out = nt)
  } else {
    if (nt != nlongitude || nt != nlatitude || nt != nuseRefraction) {
      stop("lengths of longitude, latitude and useRefraction must equal 1 or the length of time")
    }
  }
  # Ensure that the timezone is UTC. Note that Sys.Date() gives a NULL tzone.
  tzone <- attr(as.POSIXct(t[1]), "tzone")
  if (is.null(tzone) || "UTC" != tzone) {
    attributes(t)$tzone <- "UTC"
  }
  ok <- is.finite(t)
  if (!all(ok)) {
    warning("removing ", sum(!ok), " data, for which time is not finite")
    t <- t[ok]
    latitude <- latitude[ok]
    longitude <- longitude[ok]
    useRefraction <- useRefraction[ok]
  }
  # the code below is derived from fortran code, downloaded 2009-11-1 from
  # ftp://climate1.gsfc.nasa.gov/wiscombe/Solar_Rad/SunAngles/sunae.f
  t <- as.POSIXlt(t) # use this so we can work on hours, etc
  if ("UTC" != attr(as.POSIXct(t[1]), "tzone")) {
    stop("t must be in UTC")
  }
  year <- t$year + 1900
  if (any(year < 1950) || any(year > 2050)) {
    warning("year=", year[year < 1950 | year > 2050][1], " (and possibly others) is outside the acceptable range of 1950-2050")
  }
  day <- t$yday + 1
  if (any(day < 1) || any(day > 366)) {
    stop("day is not in range 1 to 366")
  }
  hour <- t$hour + t$min / 60 + t$sec / 3600
  if (any(hour < -13) || any(hour > 36)) {
    stop("hour outside range -13 to 36")
  }
  if (any(latitude < -90)) {
    warning("latitude(s) trimmed to range -90 to 90")
    latitude[latitude < -90] <- -90
  }
  if (any(latitude > 90)) {
    warning("latitude(s) trimmed to range -90 to 90")
    latitude[latitude > 90] <- 90
  }
  if (any(longitude < -180)) {
    warning("longitude(s) trimmed to range -180 to 180")
    longitude[longitude < -180] <- -180
  }
  if (any(longitude > 180)) {
    warning("longitude(s) trimmed to range -180 to 180")
    longitude[longitude > 180] <- 180
  }
  delta <- year - 1949
  leap <- delta %/% 4
  # IXME: using fortran-style int and mod here; must check for leap-year cases
  jd <- 32916.5 + (delta * 365 + leap + day) + hour / 24
  jd <- jd + ifelse(0 == (year %% 100) & 0 != (year %% 400), 1, 0)
  time <- jd - 51545
  mnlong <- 280.460 + 0.9856474 * time
  mnlong <- mnlong %% 360
  mnlong <- mnlong + ifelse(mnlong < 0, 360, 0)
  mnanom <- 357.528 + 0.9856003 * time
  mnanom <- mnanom %% 360
  mnanom <- mnanom + ifelse(mnanom < 0, 360, 0)
  rpd <- pi / 180
  mnanom <- mnanom * rpd
  eclong <- mnlong + 1.915 * sin(mnanom) + 0.020 * sin(2 * mnanom)
  eclong <- eclong %% 360
  eclong <- eclong + ifelse(eclong < 0, 360, 0)
  oblqec <- 23.439 - 0.0000004 * time
  eclong <- eclong * rpd
  oblqec <- oblqec * rpd
  num <- cos(oblqec) * sin(eclong)
  den <- cos(eclong)
  ra <- atan(num / den)
  ra <- ra + ifelse(den < 0, pi, ifelse(num < 0, 2 * pi, 0))
  dec <- asin(sin(oblqec) * sin(eclong))
  gmst <- 6.697375 + 0.0657098242 * time + hour
  gmst <- gmst %% 24
  gmst <- gmst + ifelse(gmst < 0, 24, 0)
  lmst <- gmst + longitude / 15
  lmst <- lmst %% 24
  lmst <- lmst + ifelse(lmst < 0, 24, 0)
  lmst <- lmst * 15 * rpd
  ha <- lmst - ra
  ha <- ha + ifelse(ha < (-pi), 2 * pi, 0)
  ha <- ha - ifelse(ha > pi, 2 * pi, 0)
  el <- asin(sin(dec) * sin(latitude * rpd) + cos(dec) * cos(latitude * rpd) * cos(ha))
  # pin the arg to range -1 to 1 (issue 1004)
  sinAz <- -cos(dec) * sin(ha) / cos(el)
  az <- ifelse(sinAz < (-1), -pi / 2,
               ifelse(sinAz > 1, pi / 2,
                      asin(sinAz)
               )
  )
  az <- ifelse(sin(dec) - sin(el) * sin(latitude * rpd) > 0,
               ifelse(sin(az) < 0, az + 2 * pi, az),
               pi - az
  )
  el <- el / rpd
  az <- az / rpd
  el <- el + ifelse(useRefraction,
                    ifelse(el >= 19.225,
                           0.00452 * 3.51823 / tan(el * rpd),
                           ifelse(el > (-0.766) & el < 19.225,
                                  3.51823 * (0.1594 + el * (0.0196 + 0.00002 * el)) / (1 + el * (0.505 + 0.0845 * el)),
                                  0
                           )
                    ),
                    0
  )
  soldst <- 1.00014 - 0.01671 * cos(mnanom) - 0.00014 * cos(2 * mnanom)
  soldia <- 0.5332 / soldst
  sunDRA <- oce::sunDeclinationRightAscension(t, apparent = FALSE)
  list(
    time = t, azimuth = az, altitude = el, diameter = soldia, distance = soldst,
    declination = sunDRA$declination, rightAscension = sunDRA$rightAscension
  )
}

moonAngle <- function(t, longitude = 0, latitude = 0, useRefraction = TRUE) {
  if (missing(t)) {
    stop("must provide 't'")
  }
  if (is.character(t)) {
    t <- as.POSIXct(t, tz = "UTC")
  } else if (inherits(t, "Date")) {
    t <- as.POSIXct(t)
  }
  if (!inherits(t, "POSIXt")) {
    if (is.numeric(t)) {
      tref <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC") # arbitrary
      t <- t - as.numeric(tref) + tref
    } else {
      stop("t must be POSIXt or a number corresponding to POSIXt (in UTC)")
    }
  }
  # Ensure that the timezone is UTC. Note that Sys.Date() gives a NULL tzone.
  tzone <- attr(as.POSIXct(t[1]), "tzone")
  if (is.null(tzone) || "UTC" != tzone) {
    attributes(t)$tzone <- "UTC"
  }
  RPD <- atan2(1.0, 1.0) / 45.0 # radians per degree
  # In this code, the symbol names follow Meeus (1982) chapter 30, with e.g. "p"
  # used to indicate primes, e.g. Lp stands for L' in Meeus' notation.
  # Also, T2 and T3 are powers on T.
  T <- oce::julianCenturyAnomaly(julianDay(t)) # nolint
  T2 <- T * T # nolint
  T3 <- T * T2 # nolint
  # Step 1 (top of Meuus page 148, chapter 30): mean quantities
  # moon mean longitude
  Lp <- 270.434164 + 481267.8831 * T - 0.001133 * T2 + 0.0000019 * T3 # nolint
  # sun mean amomaly
  M <- 358.475833 + 35999.0498 * T - 0.000150 * T2 - 0.0000033 * T3 # nolint
  # moon mean amomaly
  Mp <- 296.104608 + 477198.8491 * T + 0.009192 * T2 + 0.0000144 * T3 # nolint
  # moon mean elongation
  D <- 350.737486 + 445267.1142 * T - 0.001436 * T2 + 0.0000019 * T3 # nolint
  # moon distance from ascending node
  F <- 11.250889 + 483202.0251 * T - 0.003211 * T2 - 0.0000003 * T3 # nolint
  # longitude of moon ascending node
  Omega <- 259.183275 - 1934.1420 * T + 0.002078 * T2 + 0.0000022 * T3 # nolint
  # Step 2 (to bottom of p 148, chapter 30): add periodic variations ("additive terms")
  # note that 'tmp' is redefined every few lines
  tmp <- sin(RPD * (51.2 + 20.2 * T)) # nolint
  Lp <- Lp + 0.000233 * tmp
  M <- M - 0.001778 * tmp
  Mp <- Mp + 0.000817 * tmp
  D <- D + 0.002011 * tmp
  tmp <- 0.003964 * sin(RPD * (346.560 + 132.870 * T - 0.0091731 * T2)) # nolint
  Lp <- Lp + tmp
  Mp <- Mp + tmp
  D <- D + tmp
  F <- F + tmp # nolint
  tmp <- sin(RPD * Omega)
  Lp <- Lp + 0.001964 * tmp
  Mp <- Mp + 0.002541 * tmp
  D <- D + 0.001964 * tmp
  F <- F - 0.024691 * tmp # nolint
  F <- F - 0.004328 * sin(RPD * (Omega + 275.05 - 2.30 * T)) # nolint
  # Step 3: Meeus p 149
  e <- 1 - 0.002495 * T - 0.00000752 * T2 # nolint
  e2 <- e * e
  lambda <- Lp + # nolint
    (6.288750 * sin(RPD * (Mp))) + # nolint
    (1.274018 * sin(RPD * (2 * D - Mp))) + # nolint
    (0.658309 * sin(RPD * (2 * D))) + # nolint
    (0.213616 * sin(RPD * (2 * Mp))) + # nolint
    (e * -0.185596 * sin(RPD * (M))) + # nolint
    (-0.114336 * sin(RPD * (2 * F))) + # nolint
    (0.058793 * sin(RPD * (2 * D - 2 * Mp))) + # nolint
    (e * 0.057212 * sin(RPD * (2 * D - M - Mp))) + # nolint
    (0.053320 * sin(RPD * (2 * D + Mp))) + # nolint
    (e * 0.045874 * sin(RPD * (2 * D - M))) + # nolint
    (e * 0.041024 * sin(RPD * (Mp - M))) + # nolint
    (-0.034718 * sin(RPD * (D))) + # nolint
    (e * -0.030465 * sin(RPD * (M + Mp))) + # nolint
    (0.015326 * sin(RPD * (2 * D - 2 * F))) + # nolint
    (-0.012528 * sin(RPD * (2 * F + Mp))) + # nolint
    (-0.010980 * sin(RPD * (2 * F - Mp))) + # nolint
    (0.010674 * sin(RPD * (4 * D - Mp))) + # nolint
    (0.010034 * sin(RPD * (3 * M))) + # nolint
    (0.008548 * sin(RPD * (4 * D - 2 * Mp))) + # nolint
    (e * -0.007910 * sin(RPD * (M - Mp + 2 * D))) + # nolint
    (e * -0.006783 * sin(RPD * (2 * D + M))) + # nolint
    (0.005162 * sin(RPD * (Mp - D))) + # nolint
    (e * 0.005000 * sin(RPD * (M + D))) + # nolint
    (e * 0.004049 * sin(RPD * (Mp - M + 2 * D))) + # nolint
    (0.003996 * sin(RPD * (2 * Mp + 2 * D))) + # nolint
    (0.003862 * sin(RPD * (4 * D))) + # nolint
    (0.003665 * sin(RPD * (2 * D - 3 * Mp))) + # nolint
    (e * 0.002696 * sin(RPD * (2 * Mp - M))) + # nolint
    (0.002602 * sin(RPD * (Mp - 2 * F - 2 * D))) + # nolint
    (e * 0.002396 * sin(RPD * (2 * D - M - 2 * Mp))) + # nolint
    (-0.002349 * sin(RPD * (Mp + D))) + # nolint
    (e2 * 0.002249 * sin(RPD * (2 * D - 2 * M))) + # nolint
    (e * -0.002125 * sin(RPD * (2 * Mp + M))) + # nolint
    (e2 * -0.002079 * sin(RPD * (2 * M))) + # nolint
    (e2 * 0.002059 * sin(RPD * (2 * D - Mp - 2 * M))) + # nolint
    (-0.001773 * sin(RPD * (Mp + 2 * D - 2 * F))) + # nolint
    (-0.001595 * sin(RPD * (2 * F + 2 * D))) + # nolint
    (e * 0.001220 * sin(RPD * (4 * D - M - Mp))) + # nolint
    (-0.001110 * sin(RPD * (2 * Mp + 2 * F))) + # nolint
    (0.000892 * sin(RPD * (Mp - 3 * D))) + # nolint
    (e * -0.000811 * sin(RPD * (M + Mp + 2 * D))) + # nolint
    (e * 0.000761 * sin(RPD * (4 * D - M - 2 * Mp))) + # nolint
    (e2 * 0.000717 * sin(RPD * (Mp - 2 * M))) + # nolint
    (e2 * 0.000704 * sin(RPD * (Mp - 2 * M - 2 * D))) + # nolint
    (e * 0.000693 * sin(RPD * (M - 2 * Mp + 2 * D))) + # nolint
    (e * 0.000598 * sin(RPD * (2 * D - M - 2 * F))) + # nolint
    (0.000550 * sin(RPD * (Mp + 4 * D))) + # nolint
    (0.000538 * sin(RPD * (4 * Mp))) + # nolint
    (e * 0.000521 * sin(RPD * (4 * D - M))) + # nolint
    (0.000486 * sin(RPD * (2 * M - D))) # nolint
  lambda <- lambda %% 360
  B <- 0 + # nolint
    (5.128189 * sin(RPD * (F))) + # nolint
    (0.280606 * sin(RPD * (Mp + F))) + # nolint
    (0.277693 * sin(RPD * (Mp - F))) + # nolint
    (0.173238 * sin(RPD * (2 * D - F))) + # nolint
    (0.055413 * sin(RPD * (2 * D + F - Mp))) + # nolint
    (0.046272 * sin(RPD * (2 * D - F - Mp))) + # nolint
    (0.032573 * sin(RPD * (2 * D + F))) + # nolint
    (0.017198 * sin(RPD * (2 * Mp + F))) + # nolint
    (0.009267 * sin(RPD * (2 * D + Mp - F))) + # nolint
    (0.008823 * sin(RPD * (2 * Mp - F))) + # nolint
    (0.008247 * sin(RPD * (2 * D - M - F))) + # nolint
    (0.004323 * sin(RPD * (2 * D - F - 2 * Mp))) + # nolint
    (0.004200 * sin(RPD * (2 * D + F + Mp))) + # nolint
    (e * 0.003372 * sin(RPD * (F - M - 2 * D))) + # nolint
    (e * 0.002472 * sin(RPD * (2 * D + F - M - Mp))) + # nolint
    (e * 0.002222 * sin(RPD * (2 * D + F - M))) + # nolint
    (e * 0.002072 * sin(RPD * (2 * D - F - M - Mp))) + # nolint
    (e * 0.001877 * sin(RPD * (F - M + Mp))) + # nolint
    (0.001828 * sin(RPD * (4 * D - F - Mp))) + # nolint
    (e * -0.001803 * sin(RPD * (F + M))) + # nolint
    (-0.001750 * sin(RPD * (3 * F))) + # nolint
    (e * 0.001570 * sin(RPD * (Mp - M - F))) + # nolint
    (-0.001487 * sin(RPD * (F + D))) + # nolint
    (e * -0.001481 * sin(RPD * (F + M + Mp))) + # nolint
    (e * 0.001417 * sin(RPD * (F - M - Mp))) + # nolint
    (e * 0.001350 * sin(RPD * (F - M))) + # nolint
    (0.001330 * sin(RPD * (F - D))) + # nolint
    (0.001106 * sin(RPD * (F + 3 * Mp))) + # nolint
    (0.001020 * sin(RPD * (4 * D - F))) + # nolint
    (0.000833 * sin(RPD * (F + 4 * D - Mp))) + # nolint
    (0.000781 * sin(RPD * (Mp - 3 * F))) + # nolint
    (0.000670 * sin(RPD * (F + 4 * D - 2 * Mp))) + # nolint
    (0.000606 * sin(RPD * (2 * D - 3 * F))) + # nolint
    (0.000597 * sin(RPD * (2 * D + 2 * Mp - F))) + # nolint
    (e * 0.000492 * sin(RPD * (2 * D + Mp - M - F))) + # nolint
    (0.000450 * sin(RPD * (2 * Mp - F - 2 * D))) + # nolint
    (0.000439 * sin(RPD * (3 * Mp - F))) + # nolint
    (0.000423 * sin(RPD * (F + 2 * D + 2 * Mp))) + # nolint
    (0.000422 * sin(RPD * (2 * D - F - 3 * Mp))) + # nolint
    (e * -0.000367 * sin(RPD * (F + F + 2 * D - Mp))) + # nolint
    (e * -0.000353 * sin(RPD * (M + F + 2 * D))) + # nolint
    (0.000331 * sin(RPD * (F + 4 * D))) + # nolint
    (e * 0.000317 * sin(RPD * (2 * D + F - M + Mp))) + # nolint
    (e2 * 0.000306 * sin(RPD * (2 * D - 2 * M - F))) + # nolint
    (-0.000283 * sin(RPD * (Mp + 3 * F))) # nolint
  omega1 <- 0.0004664 * cos(RPD * Omega)
  omega2 <- 0.0000754 * cos(RPD * (Omega + 275.05 - 2.30 * T)) # nolint
  beta <- B * (1 - omega1 - omega2)
  pi <- 0.950724 + # nolint
    (0.051818 * cos(RPD * (Mp))) + # nolint
    (0.009531 * cos(RPD * (2 * D - Mp))) + # nolint
    (0.007843 * cos(RPD * (2 * D))) + # nolint
    (0.002824 * cos(RPD * (2 * Mp))) + # nolint
    (0.000857 * cos(RPD * (2 * D + Mp))) + # nolint
    (e * 0.000533 * cos(RPD * (2 * D - M))) + # nolint
    (e * 0.000401 * cos(RPD * (2 * D - M - Mp))) + # nolint
    (e * 0.000320 * cos(RPD * (Mp - M))) + # nolint
    (-0.000271 * cos(RPD * (D))) + # nolint
    (e * -0.000264 * cos(RPD * (M + Mp))) + # nolint
    (-0.000198 * cos(RPD * (2 * F - Mp))) + # nolint
    (0.000173 * cos(RPD * (3 * Mp))) + # nolint
    (0.000167 * cos(RPD * (4 * D - Mp))) + # nolint
    (e * -0.000111 * cos(RPD * (M))) + # nolint
    (0.000103 * cos(RPD * (4 * D - 2 * Mp))) + # nolint
    (-0.000084 * cos(RPD * (2 * Mp - 2 * D))) + # nolint
    (e * -0.000083 * cos(RPD * (2 * D + M))) + # nolint
    (0.000079 * cos(RPD * (2 * D + 2 * Mp))) + # nolint
    (0.000072 * cos(RPD * (4 * D))) + # nolint
    (e * 0.000064 * cos(RPD * (2 * D - M + Mp))) + # nolint
    (e * -0.000063 * cos(RPD * (2 * D + M - Mp))) + # nolint
    (e * 0.000041 * cos(RPD * (M + D))) + # nolint
    (e * 0.000035 * cos(RPD * (2 * Mp - M))) + # nolint
    (-0.000033 * cos(RPD * (3 * Mp - 2 * D))) + # nolint
    (-0.000030 * cos(RPD * (Mp + D))) + # nolint
    (-0.000029 * cos(RPD * (2 * F - 2 * D))) + # nolint
    (e * -0.000029 * cos(RPD * (2 * Mp + M))) + # nolint
    (e2 * 0.000026 * cos(RPD * (2 * D - 2 * M))) + # nolint
    (-0.000023 * cos(RPD * (2 * F - 2 * D + Mp))) + # nolint
    (e * 0.000019 * cos(RPD * (4 * D - M - Mp))) # nolint
  # For coordinate conversions, need epsilon (obliquity of the ecliptic)
  # as defined in Meuus eq 18.4, page 81.
  epsilon <- 23.452294 - 0.0130125 * T - 0.00000164 * T2 + 0.000000503 * T3 # nolint
  ec <- eclipticalToEquatorial(lambda, beta, epsilon)
  # .lh <- equatorialToLocalHorizontal(ec$rightAscension, ec$declination, t, latitude, longitude)
  lh <- equatorialToLocalHorizontal(
    rightAscension = ec$rightAscension,
    declination = ec$declination,
    t = t,
    longitude = longitude,
    latitude = latitude
  )
  # Illuminated fraction, reference 1 chapter 31 (second, approximate, formula)
  D <- D %% 360 # need this; could have done it earlier, actually
  illfr <- 180 - D - 6.289 * sin(RPD * Mp) +
    2.100 * sin(RPD * M) -
    1.274 * sin(RPD * (2 * D - Mp)) -
    0.658 * sin(RPD * 2 * D) -
    0.2114 * sin(RPD * 2 * Mp) -
    0.112 * sin(RPD * D)
  illuminatedFraction <- (1 + cos(RPD * illfr)) / 2
  # Meeus (1982) eq 32.3 page 160
  phase <- T * 1236.85 # nolint
  
  # The 180 in azimuth converts from astron convention with azimuth=westward
  # angle from South, to eastward from North.
  res <- list(
    time = t,
    azimuth = lh$azimuth + 180,
    altitude = lh$altitude,
    rightAscension = ec$rightAscension, declination = ec$declination,
    lambda = lambda %% 360, beta = beta,
    diameter = pi, distance = 6378.14 / sin(RPD * pi),
    illuminatedFraction = illuminatedFraction,
    phase = phase
  )
  res
}
