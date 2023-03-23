#' Pivot Ornitela accelerometry data to long format (xyzxyzxyz...)
#'
#' @description This function takes an data frame made from a .csv downloaded
#'   directly from the Ornitela website, and pivots it into the format
#'   "xyzxyzxyz...", which is used for analyses in the package \code{rabc} on
#'   Github.
#'
#' @param data A data frame of an Ornitela .csv for acc data ("SENSORS").
#' @param burstMinApart The frequency of ACC bursts in minutes.
#' @param pointsPerBurst The number of points in a burst, a.k.a. hertz * seconds.
#' @param removeIncompleteBursts Should incomplete bursts (bursts with fewer
#'   than \code{pointsPerBurst}) be removed from the final data frame?
#'
#' @return A data frame.
#' @export
#'
#' @examples
#' NA
acc_pivot_wider <- function(data, burstMinApart = 6, pointsPerBurst = 30,
                      removeIncompleteBursts = TRUE) {
  # try to update this to auto-detect burst statistics,
  # rather than force them to be entered as function arguments
  data <- data %>%
    filter(., .$datatype == "SENSORS") %>%
    select(device_id, UTC_datetime, acc_x, acc_y, acc_z)
  data$UTC_datetime <- as.POSIXct(data$UTC_datetime,
                                  tryFormats = c("%m/%d/%y %H:%M:%S",
                                                 "%m/%d/%y %H:%M",
                                                 "%m/%d/%Y %H:%M",
                                                 "%m/%d/%Y %H:%M:%S"), tz = "GMT")
  uniqueBursts <- unique(data$UTC_datetime)
  bursts<- c(uniqueBursts[1])
  for (i in 2:length(uniqueBursts)){
    if (difftime(uniqueBursts[i],
                 uniqueBursts[i-1], units = "mins") > (burstMinApart/2)){
      burst <- c(uniqueBursts[i])
      bursts <- append(bursts,burst)
    }
  }
  result <- data.frame(matrix(NA,
                              nrow = length(bursts),
                              ncol = (3*pointsPerBurst)+2))
  colnames(result)[1] <- "device_id"
  colnames(result)[2] <- "UTC_datetime"
  colnames(result)[3:ncol(result)] <- c("x","y","z")
  result[,1] <- data$device_id[1]
  result[,2] <- bursts
  vec <- c(data$acc_x[1],data$acc_y[1],data$acc_z[1])
  burstPosition <- 1
  for (i in 2:length(data$acc_x)){
    if (difftime(data$UTC_datetime[i],
                 data$UTC_datetime[i-1],
                 units = "mins") < (burstMinApart/2)){
      vec.new <- c(data$acc_x[i],data$acc_y[i],data$acc_z[i])
      vec <- append(vec,vec.new)
      if (i == length(data$acc_x)){
        len <- length(vec)
        result[burstPosition, 3:(len+2)] <- vec
      }
    } else {
      len <- length(vec)
      result[burstPosition, 3:(len+2)] <- vec
      burstPosition <- burstPosition + 1
      vec <- c(data$acc_x[i],data$acc_y[i],data$acc_z[i])
    }
  }
  if (removeIncompleteBursts == TRUE){
    result <- result[complete.cases(result),]
  }
  return(result)
}
