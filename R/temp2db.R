#' tmp2db
#'
#' @param year 
#' @param haul 
#' @param hachi 
#' @param sets 
#'
#' @return
#' @export
#'
#' @examples
temp2db <- function(year = year, haul = haul, hachi = hachi, sets = sets) {
  # Setup to record 0 to 820 meters, by 1-m increment
  Depth <- seq(0, 1200, 1)
  
  # Create data frames to store data
  StDat <- data.frame() # For 1-m increment temperature profile
  cleanBot <- data.frame() # 
  cleanProf <- data.frame()
  
  for (ab in c("", "b")) {
    for (sta in c(1:535)) {
      if (!file.exists(paste(yr, "/sta", sta, ab, ".asc", sep = ""))) {
        next
      }
      
      # # Pull haul data for given year/station
      hal <- haul %>% dplyr::filter(station == sta)
      
      # Read in data from tdr
      tdr <- read.table(paste(yr, "/sta", sta, ab, ".asc", sep = ""), skip = 49, header = FALSE, sep = ",")
      colnames(tdr) <- c("Temperature", "Depth", "Date", "Time")
      
      tdr$Diff <- c(0, diff(tdr$Depth))
      tdr$Time <- as.POSIXct(strptime(tdr$Time, format = "%H:%M:%S"))
      tdr$Julian <- as.numeric(format(as.Date(tdr$Date, format = "%d %b %Y"), "%j"))
      
      # Bottom depth and temp from TDR, stationary
      bot <- as.data.frame(tdr[abs(tdr$Diff) <= 0.03 & tdr$Depth > (max(tdr$Depth) - 30), ])
      bot <- bot[abs(bot$Depth - mean(bot$Depth)) <= 5, ]
      
      Year <- yr
      Julian <- round(mean(bot$Julian), 0)
      Station <- sta
      TDRDepth <- round(mean(bot$Depth), 1)
      minDepth <- round(min(bot$Depth), 1)
      maxDepth <- round(max(bot$Depth), 1)
      MeanTemp <- round(mean(bot$Temperature), 2)
      MinTemp <- round(min(bot$Temperature), 2)
      MaxTemp <- round(max(bot$Temperature), 2)
      NumScans <- nrow(bot)
      
      # Downcast
      dodat <- tdr[tdr$Diff > 1.2 & tdr$Depth > 0, ]
      dodat <- rbind(dodat, bot[bot$Time == min(bot$Time), ])
      dodat$TimeDiff <- dodat$Time - mean(dodat$Time)
      units(dodat$TimeDiff) <- "mins"
      dodat <- dodat[abs(dodat$TimeDiff) < 60, ]
      
      Temp1m <- oce::oce.approx(dodat$Depth, dodat$Temperature, Depth, "rr")
      newD <- data.frame(cbind(Depth,Temp1m))
      newD$Station <- sta
      newD$Year <- yr
      newD$Set <- ifelse (ab == "b", "Two", "One")
      newDat <- na.omit(newD)
      
      StDat <- rbind(StDat, newDat)
      
      # Upcast
      updat <- tdr[tdr$Diff < -1.2 & tdr$Depth > 0, ]
      updat$TimeDiff <- updat$Time - mean(updat$Time)
      units(updat$TimeDiff) <- "mins"
      updat <- updat[abs(updat$TimeDiff) < 40, ]
      
      # Read in data for each hachi
      hach <- hachi[hachi$station == sta & !is.na(hachi$time), ]
      
      # If there is an upcast, get the time that the TDR hits the surface (ie, when that hachi is recorded by observer at the rail)
      TDRtime <- ifelse (nrow(updat) == 0, NA, as.numeric(paste0(str_pad(hour(updat[nrow(updat), 4]), 
                                                                         2, pad="0"), str_pad(minute(updat[nrow(updat), 4]), 2, pad="0"))))
      
      # Using the TDR estimated time at the rail, estimate the hachi's time using a 10 minute window
      HachiTimeEst <- ceiling(ifelse (is.na(TDRtime) || TDRtime < min(hach$time) || 
                                        TDRtime > max(hach$time), NA, median(hach[!is.na(hach$time) & 
                                                                                    hach$time < as.integer(paste0(str_pad(hour(updat[nrow(updat), 4]), 2, pad = "0"),
                                                                                                                  str_pad(minute(updat[nrow(updat), 4]), 2, pad = "0"))) + 5 & 
                                                                                    hach$time > as.integer(paste0(str_pad(hour(updat[nrow(updat), 4]), 2, pad = "0"),
                                                                                                                  str_pad(minute(updat[nrow(updat), 4]), 2, pad = "0"))) - 5, "hachi"])))
      
      # From the time estimated above, assign an estimated hachi number, and when that is NA, assign to skate 60, the presumed skate number
      HachiNumEst <- ifelse (is.na(HachiTimeEst), 60, HachiTimeEst)
      HachiDefTime <- ifelse (min(hach$hachi) > 60, NA, hach[hach$hachi == 60, "time"]) # Time at hachi 60, the presumed skate
      DepthExp <- ifelse (min(hach$hachi) > 60, NA, hach[hach$hachi == 60, "intrpdep"]) # Depth at presumed skate
      
      # Pull the observer assigned depth (an interpolated value from the skates they actually record) at the skate believed to have the TDR
      ObsDepth <- hach[hach$hachi == HachiNumEst, "intrpdep"]
      
      # Estimated hachi number from time estimate using TDR timestamp
      hachTDR <- data.frame(cbind(HachiNumEst, TDRDepth))
      names(hachTDR) <- c("hachi", "depth")
      
      # Select haul, in case of 2-haul station, just morning, one haul station it is the only haul at the station, so could be afternoon with hachi numbers not starting at 1
      haulOne <- hal[hal$haul == min(hach$haul), ]
      
      # Get depth bounds from haul sheet and put on either end of the set, not perfect, but good for guidance
      Zx1 <- as.data.frame(c(min(hach[hach$haul == haulOne[, "haul"], "hachi"]) - 0.5, max(hach[hach$haul == haulOne[, "haul"], "hachi"]) + 0.5))
      Zy1 <- as.data.frame(t(haulOne[, c("starting_depth", "ending_depth")]))
      one <- cbind(Zx1, Zy1)
      names(one) <- c("hachi", "depth")
      
      chx <- rbind(one, hachTDR)
      
      # Used for estimating location between start and end of set as percent
      row <- round(100 * ((hachTDR$hachi - min(hach[hach$haul == haulOne[, "haul"], "hachi"]) - 1) / (max(hach[hach$haul == haulOne[, "haul"], "hachi"]) - min(hach[hach$haul == haulOne[, "haul"], "hachi"]) - 1)), 0)
      
      # Get Lat and Lon from where TDR is guessed to be, first in lat/lon then in UTM for plots in km assuming straight line between points
      Lat <- haulOne$start.lat + row  / 100 * (haulOne$end.lat - haulOne$start.lat)
      Lon <- haulOne$start.lon + row / 100 * (haulOne$end.lon - haulOne$start.lon)
      
      # Get UTM coordinates of start and end of set, used for visualizing errors
      pt <- data.frame(cbind(Lon, Lat))
      UTMzn <- floor((Lon + 180) / 6) + 1
      setp <- SpatialPoints(pt, proj4string = CRS("+proj=longlat"))
      
      ptsOnest <- data.frame(spTransform(SpatialPoints(haulOne[, c("start.lon", "start.lat")], proj4string = CRS("+proj=longlat")), 
                                         CRS(paste0("+proj=utm +zone=", UTMzn, " +datum=WGS84 +units=km"))))
      ptsOneen <- data.frame(spTransform(SpatialPoints(haulOne[, c("end.lon", "end.lat")], proj4string = CRS("+proj=longlat")),
                                         CRS(paste0("+proj=utm +zone=", UTMzn, " +datum=WGS84 +units=km"))))
      
      # If there are more than one haul at a station, such as normal two-haul station
      if (length(unique(hach$haul)) > 1) {
        haulTwo <- hal[!hal$haul==min(hach$haul),]
        
        Zx2 <- as.data.frame(c(max(hach$hachi) / 2 + 0.75, max(hach$hachi) + 0.5))
        Zy2 <- as.data.frame(t(haulTwo[, c("starting_depth", "ending_depth")]))
        two <- cbind(Zx2, Zy2)
        names(two) <- c("hachi", "depth")
        
        chx <- rbind(one, two, hachTDR)
        
        ptsTwost <- data.frame(spTransform(SpatialPoints(haulTwo[, c("start.lon", "start.lat")], proj4string = CRS("+proj=longlat")), 
                                           CRS(paste0("+proj=utm +zone=", UTMzn, " +datum=WGS84 +units=km"))))
        ptsTwoen <- data.frame(spTransform(SpatialPoints(haulTwo[, c("end.lon", "end.lat")], proj4string = CRS("+proj=longlat")),
                                           CRS(paste0("+proj=utm +zone=", UTMzn, " +datum=WGS84 +units=km"))))
        
        if (hachTDR$hachi > 90) {
          row <- round(100 * ((hachTDR$hachi - min(hach[hach$haul == haulTwo[, "haul"], "hachi"]) - 1) / (max(hach[hach$haul == haulTwo[, "haul"], "hachi"]) - min(hach[hach$haul == haulTwo[, "haul"], "hachi"]) - 1)), 0)
          
          Lat <- haulTwo$start.lat + row / 100 * (haulTwo$end.lat - haulTwo$start.lat)
          Lon <- haulTwo$start.lon + row / 100 * (haulTwo$end.lon - haulTwo$start.lon)
        }
      }
      
      UTMp <- spTransform(setp, CRS(paste0("+proj=utm +zone=", UTMzn, " +datum=WGS84 +units=km")))
      UTMdf <- data.frame(UTMp)
      Lat2m <- UTMdf[, 2]
      Lon2m <- UTMdf[, 1]
      
      # Converting to UTM, and adjusting track line to fall between actual start/end but with trajectory of master, to capture curvature
      # Not equal number of points or spacing, so create 100 segments using the points provided, then row/100 puts it at the estimated TDR along track
      # Not all stations have a reference, so they are left with their assumed straight line track
      
      if (paste0(sta, "A") %in% unique(sets$StnSet))  {
        setOnep <- SpatialPoints(sets[sets$StnSet == paste0(sta, "A"), c(3, 2)], proj4string = CRS("+proj=longlat"))
        
        UTMzn <- floor((mean(setOnep$Longitude) + 180) / 6) + 1
        setUTM <- spTransform(setOnep, CRS(paste0("+proj=utm +zone=", UTMzn, " +datum=WGS84 +units=km")))
        setOne <- Lines(list(Line(setUTM)), ID = paste0(sta, "A"))
        
        testOne <- spsample(setOne, n = 101, type = "regular")
        testOnedf <- data.frame(testOne)
        testOnedf$per <- c(1:101)
        
        ptsOnest <- data.frame(spTransform(SpatialPoints(haulOne[, c("start.lon", "start.lat")], proj4string = CRS("+proj=longlat")), 
                                           CRS(paste0("+proj=utm +zone=", UTMzn, " +datum=WGS84 +units=km"))))
        ptsOneen <- data.frame(spTransform(SpatialPoints(haulOne[, c("end.lon", "end.lat")], proj4string = CRS("+proj=longlat")),
                                           CRS(paste0("+proj=utm +zone=", UTMzn, " +datum=WGS84 +units=km"))))
        
        Xioff <- as.numeric(testOnedf$Longitude[1] - as.data.frame(ptsOnest$start.lon[1]))
        Yioff <- as.numeric(testOnedf$Latitude[1] - as.data.frame(ptsOnest$start.lat[1]))
        Xfoff <- as.numeric(testOnedf$Longitude[101] - as.data.frame(ptsOneen$end.lon[1]))
        Yfoff <- as.numeric(testOnedf$Latitude[101] - as.data.frame(ptsOneen$end.lat[1]))
        
        testOnedf$long2 <- testOnedf$Longitude - (((101 - testOnedf$per) / 100) * Xioff) - (((testOnedf$per - 1) / 100) * Xfoff)
        testOnedf$lat2 <- testOnedf$Latitude - (((101 - testOnedf$per) / 100) * Yioff) - (((testOnedf$per - 1) / 100) * Yfoff)
        
        Lat2m <- testOnedf[row, 5]
        Lon2m <- testOnedf[row, 4] 
        
        pt <- data.frame(cbind(Lon2m, Lat2m))
        setp <- SpatialPoints(pt, proj4string = CRS(paste0("+proj=utm +zone=", UTMzn, " +datum=WGS84 +units=km")))
        
        LLp <- spTransform(setp, CRS(paste0("+proj=longlat")))
        LLdf <- data.frame(LLp)
        Lat <- LLdf[, 2]
        Lon <- LLdf[, 1]
        
        if (length(unique(hach$Haul)) > 1) {
          setTwop <- SpatialPoints(sets[sets$StnSet == paste0(sta, "B"), c(3, 2)], proj4string = CRS("+proj=longlat"))
          
          UTMzn2 <- floor((mean(setTwop$Longitude) + 180) / 6) + 1
          setUTM2 <- spTransform(setTwop, CRS(paste0("+proj=utm +zone=", UTMzn2, " +datum=WGS84 +units=km")))
          setTwo <- Lines(list(Line(setUTM2)), ID = paste0(sta, "B"))
          
          testTwo <- spsample(setTwo, n = 101, type = "regular")
          testTwodf <- data.frame(testTwo)
          testTwodf$per <- c(1:101)
          
          ptsTwost <- data.frame(spTransform(SpatialPoints(haulTwo[, c("start.lon", "start.lat")], proj4string = CRS("+proj=longlat")), 
                                             CRS(paste0("+proj=utm +zone=", UTMzn, " +datum=WGS84 +units=km"))))
          ptsTwoen <- data.frame(spTransform(SpatialPoints(haulTwo[, c("end.lon", "end.lat")], proj4string = CRS("+proj=longlat")),
                                             CRS(paste0("+proj=utm +zone=", UTMzn, " +datum=WGS84 +units=km"))))
          
          Xioff2 <- as.numeric(testTwodf$Longitude[1] - as.data.frame(ptsTwost$start.lon[1]))
          Yioff2 <- as.numeric(testTwodf$Latitude[1] - as.data.frame(ptsTwost$start.lat[1]))
          Xfoff2 <- as.numeric(testTwodf$Longitude[101] - as.data.frame(ptsTwoen$end.lon[1]))
          Yfoff2 <- as.numeric(testTwodf$Latitude[101] - as.data.frame(ptsTwoen$end.lat[1]))
          
          testTwodf$long2 <- testTwodf$Longitude - (((101 - testTwodf$per) / 100) * Xioff2) - (((testTwodf$per - 1) / 100) * Xfoff2)
          testTwodf$lat2 <- testTwodf$Latitude - (((101 - testTwodf$per) / 100) * Yioff2) - (((testTwodf$per - 1) / 100) * Yfoff2)
          
          # If the TDR is only on the second haul, and the hachi numbers don't actually start at 1
          if (hachTDR$hachi > 90) {
            Lat2m <- testTwodf[row, 5]
            Lon2m <- testTwodf[row, 4] 
            
            pt <- data.frame(cbind(Lon2m, Lat2m))
            setp <- SpatialPoints(pt, proj4string = CRS(paste0("+proj=utm +zone=", UTMzn, " +datum=WGS84 +units=km")))
            
            LLp <- spTransform(setp, CRS(paste0("+proj=longlat")))
            LLdf <- data.frame(LLp)
            Lat <- LLdf[, 2]
            Lon <- LLdf[, 1]
          }
        }
      }
      
      Set <- ifelse (ab == "b", "Two", "One")
      newBot <- data.frame(Year, Station, Julian, Lat, Lon, Lat2m, Lon2m, MeanTemp, MinTemp, MaxTemp, TDRDepth, minDepth, maxDepth, TDRtime, HachiTimeEst, HachiNumEst, ObsDepth, DepthExp, HachiDefTime, Set)
      cleanBot <- rbind(cleanBot, newBot)
      
      if (nrow(dodat) > 0) {
        dodat$Year <- Year
        dodat$Station <- Station
        # dodat$Reg = Reg
        # dodat$Hab = Hab
        # dodat$Geo = Geo
        dodat$Julian <- Julian
        dodat$Lat <- Lat
        dodat$Lon <- Lon 
        dodat$Lat2m <- Lat2m
        dodat$Lon2m <- Lon2m
        dodat$maxZ <- TDRDepth
        dodat$meanBot <- MeanTemp
        dodat$minBot <- MinTemp
        dodat$maxBot <- MaxTemp
        dodat$Set <- ifelse (ab == "b", "Two", "One")
        
        cleanProf <- rbind(cleanProf, dodat)
      }
      
      if (length(unique(hach$Haul)) > 1) {
        rm(haulTwo)
        rm(testTwodf)
        rm(ptsTwost)
        rm(ptsTwoen)
      }
    }
  }
  
}
