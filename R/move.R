#' Title
#'
#' @param data A \code{Move*} object.
#' @param birds Temporary.
#' @param date Temporary.
#' @param month Integers 1-12, representing the month of the year.
#' @param year Four digit year.
#' @param DayNight Should the \code{Move*} object be subset to just daytime or
#'   just nighttime points? Choose \code{day} for day and \code{night} for night.
#' @param timeSplit Either \code{daily}, \code{weekly}, \code{monthly}, or \code{yearly}.
#'   Creates multiple \code{Move} objects at the given frequency. For example, if
#'   your data consists of two full months and you choose \code{timesplit = "daily"},
#'   the function will return a list of 60+ \code{Move} objects, one for each day.
#'
#' @return A \code{Move} object or a list of \code{Move} objects.
#' @export
#'
#' @examples
#' NA
subsetMove <- function(data, birds = NA, date = NA, month = NA, year = NA,
                       DayNight = "no", timeSplit = "none"){
  if (is.na(birds) == TRUE &&
      is.na(date) == TRUE &&
      is.na(month) == TRUE &&
      is.na(year) == TRUE &&
      DayNight == "no" &&
      timeSplit == "none"){
    message("No subset conditions were input. Returning unchanged data.")
    return(data)
  }
  if (DayNight != "day" & DayNight != "night" & DayNight != "no"){
    stop("DayNight must be 'day', 'night', or 'no'.")
  }
  if (class(data) == "MoveStack"){
    data <- move::split(data)
  }
  if (is.list(data) == TRUE){
    results <- plyr::llply(data, subsetMove, birds = birds, date = date,
                           month = month, year = year, DayNight = DayNight,
                           timeSplit = timeSplit)
    return(results)
  }
  if (is.na(birds) == FALSE){
    data <- data[names(data) %in% birds]
  }
  if (is.na(date) == FALSE){
    data <- data[as.Date(timestamps(data))==as.Date(date)]
  }
  if (is.na(month) == FALSE){
    data <- data[month(timestamps(data))==month]
  }
  if (is.na(year) == FALSE){
    data <- data[year(timestamps(data))==year]
  }
  if (DayNight == "day"){
    DN <- rep("Day", n.locs(data)-1)
    DN[maptools::solarpos(
      data[-n.locs(data)],
      timestamps(data)[-n.locs(data)])[,2] < -6 &
        maptools::solarpos(data[-1], timestamps(data)[-1])[,2] < -6] <- "Night"
    d.burst <- move::burst(x = data, f = DN)
    data <- data[d.burst@burstId == "Day"]
  }
  if (DayNight == "night"){
    DN <- rep("Day", n.locs(data)-1)
    DN[maptools::solarpos(
      data[-n.locs(data)],
      timestamps(data)[-n.locs(data)])[,2] < -6 &
        maptools::solarpos(data[-1], timestamps(data)[-1])[,2] < -6] <- "Night"
    d.burst <- move::burst(x = data, f = DN)
    data <- data[d.burst@burstId == "Night"]
  }
  if (timeSplit == "daily"){
    startDate <- date(timestamps(data))[1]
    endDate <- date(timestamps(data))[length(data)]
    dates <- seq(startDate, endDate, by=1)
    list <- list()
    for (d in dates) {
      data.day <- data[as.Date(timestamps(data), tz = tz(timestamps(data))) == d]
      list <- append(list,data.day)
    }
    names(list) <- as.character(dates)
    return(list)
  }
  if (timeSplit == "weekly"){
    startDate <- date(timestamps(data))[1]
    endDate <- date(timestamps(data))[length(data)]
    weeks <- seq(startDate, endDate, by="week")
    list <- list()
    for (i in 1:length(weeks)) {
      data.week <- data[as.Date(timestamps(data), tz = tz(timestamps(data))) %in%
                          c(as.Date(weeks[i]):as.Date(weeks[i]+6))]
      list <- append(list,data.week)
    }
    names(list) <- as.character(weeks)
    return(list)
  }
  if (timeSplit == "monthly"){
    startDate <- date(timestamps(data))[1]
    endDate <- date(timestamps(data))[length(data)]
    months <- month(seq(startDate, endDate, by="month"), label = TRUE)
    years <- year(seq(startDate, endDate, by="month"))
    list <- list()
    for (i in 1:length(months)) {
      data.month <- data[month(timestamps(data), label = TRUE) == months[i]]
      list <- append(list,data.month)
    }
    names(list) <- paste(months, years, sep = " ")
    return(list)
  }
  if (timeSplit == "yearly"){
    startDate <- date(timestamps(data))[1]
    endDate <- date(timestamps(data))[length(data)]
    years <- unique(year(seq(startDate, endDate, by="week")))
    list <- list()
    for (i in 1:length(years)) {
      data.month <- data[year(timestamps(data))==years[i]]
      list <- append(list,data.month)
    }
    names(list) <- years
    return(list)
  }
  return(data)
}
