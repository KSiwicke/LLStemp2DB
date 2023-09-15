#' tmp2db
#'
#' @param yr 
#' @param haul 
#' @param hachi 
#' @param sets
#' @param coords 
#'
#' @return
#' @export
#'
#' @examples
temp2db <- function(yr = yr, haul = haul, hachi = hachi, sets = sets, coords = coords) {
  # Setup to record 0 to 820 meters, by 1-m increment
  Depth <- seq(0, 1200, 1)
  
  # Create data frames to store data
  StDat <- data.frame() # For 1-m increment temperature profile
  cleanBot <- data.frame() # 
  cleanProf <- data.frame()
  
  for (h in min(coords$haul):max(coords$haul)) {
    coord <- coords %>% dplyr::filter(haul == h)
    sta <- coord$station
    ab <- ifelse (coord$set == "One", "", "b")
    if (!file.exists(paste(yr, "/sta", sta, ab, ".asc", sep = ""))) {
      next
    }  
    
    # # Pull haul data for given year/station
    hal <- haul %>% dplyr::filter(haul == h)
      
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
    newD$Set <- coord$set
    newDat <- na.omit(newD)
    
    StDat <- rbind(StDat, newDat)
      
    # Upcast
    updat <- tdr[tdr$Diff < -1.2 & tdr$Depth > 0, ]
    updat$TimeDiff <- updat$Time - mean(updat$Time)
    units(updat$TimeDiff) <- "mins"
    updat <- updat[abs(updat$TimeDiff) < 40, ]
      
    # # Read in data for each hachi
    hach <- hachi[hachi$haul == h, ]
    ObsDepth <- hach[hach$hachi == coord$hachi, "intrpdep"]
      
    Lat <- coord$lat_deg + coord$lat_min / 60
    Lon <- coord$lon_deg + coord$lon_min / 60

    newBot <- data.frame(Year, Station, Julian, Lat, Lon, MeanTemp, MinTemp, MaxTemp, TDRDepth, minDepth, maxDepth, ObsDepth, Set)
    cleanBot <- rbind(cleanBot, newBot)
    
    if (nrow(dodat) > 0) {
      dodat$Year <- Year
      dodat$Station <- Station
      dodat$Julian <- Julian
      dodat$Lat <- Lat
      dodat$Lon <- Lon 
      dodat$maxZ <- TDRDepth
      dodat$meanBot <- MeanTemp
      dodat$minBot <- MinTemp
      dodat$maxBot <- MaxTemp
      dodat$Set <- Set
      
      cleanProf <- rbind(cleanProf, dodat)
    }
  }
  return(list(cleanBot, cleanProf))
}