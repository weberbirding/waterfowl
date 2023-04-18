#' Subset or split a Move* object
#'
#' @description A function to simplify subsetting or splitting up a \code{Move}
#'   object. Allows for pulling out a single month and/or year, pulling out day
#'   vs. night, and/or splitting into smaller \code{Move} objects for each day,
#'   week, month, or year.
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
#'
#' @import move
subsetMove <- function(data, birds = NA, date = NA, month = NA, year = NA,
                       DayNight = "no", timeSplit = "none"){
  if (is.na(birds) &&
      is.na(date) &&
      is.na(month) &&
      is.na(year) &&
      DayNight == "no" &&
      timeSplit == "none") {
    message("No subset conditions were input. Returning unchanged data.")
    return(data)
  }
  if (DayNight != "day" & DayNight != "night" & DayNight != "no"){
    stop("DayNight must be 'day', 'night', or 'no'.")
  }
  if (is(data, "MoveStack")) {
    data <- move::split(data)
  }
  if (is.list(data)) {
    results <- plyr::llply(data, subsetMove, birds = birds, date = date,
                           month = month, year = year, DayNight = DayNight,
                           timeSplit = timeSplit)
    return(results)
  }
  if (!is.na(birds)) {
    data <- data[names(data) %in% birds]
  }
  if (!is.na(date)) {
    data <- data[as.Date(move::timestamps(data))==as.Date(date)]
  }
  if (!is.na(month)) {
    data <- data[month(move::timestamps(data))==month]
  }
  if (!is.na(year)) {
    data <- data[year(move::timestamps(data))==year]
  }
  if (DayNight == "day") {
    tod <- time_of_day(move)
    data <- data[tod == "day"]
  }
  if (DayNight == "night") {
    tod <- time_of_day(move)
    data <- data[tod == "night"]
  }
  if (timeSplit == "daily") {
    startDate <- lubridate::date(move::timestamps(data)[1])
    endDate <- lubridate::date(move::timestamps(data)[length(data)])
    dates <- seq(startDate, endDate, by = 1)
    list <- list()
    for (d in dates) {
      data.day <- data[as.Date(move::timestamps(data),
                               tz = lubridate::tz(move::timestamps(data))) == d]
      list <- append(list, data.day)
    }
    names(list) <- as.character(dates)
    return(list)
  }
  if (timeSplit == "weekly") {
    startDate <- lubridate::date(move::timestamps(data)[1])
    endDate <- lubridate::date(move::timestamps(data)[length(data)])
    weeks <- seq(startDate, endDate, by="week")
    list <- list()
    for (i in 1:length(weeks)) {
      data.week <- data[as.Date(move::timestamps(data),
                                tz = lubridate::tz(move::timestamps(data))) %in%
                          c(as.Date(weeks[i]):as.Date(weeks[i] + 6))]
      list <- append(list,data.week)
    }
    names(list) <- as.character(weeks)
    return(list)
  }
  if (timeSplit == "monthly") {
    startDate <- lubridate::date(move::timestamps(data)[1])
    endDate <- lubridate::date(move::timestamps(data)[length(data)])
    months <- lubridate::month(seq(startDate, endDate, by = "month"), label = TRUE)
    years <- lubridate::year(seq(startDate, endDate, by = "month"))
    list <- list()
    for (i in 1:length(months)) {
      data.month <- data[lubridate::month(move::timestamps(data), label = TRUE) == months[i]]
      list <- append(list,data.month)
    }
    names(list) <- paste(months, years, sep = " ")
    return(list)
  }
  if (timeSplit == "yearly") {
    startDate <- lubridate::date(move::timestamps(data)[1])
    endDate <- lubridate::date(move::timestamps(data)[length(data)])
    years <- unique(lubridate::year(seq(startDate, endDate, by = "week")))
    list <- list()
    for (i in 1:length(years)) {
      data.month <- data[lubridate::year(move::timestamps(data)) == years[i]]
      list <- append(list, data.month)
    }
    names(list) <- years
    return(list)
  }
  return(data)
}

#' Find the time of day
#'
#' @description A simple wrapper around \code{suncalc::getSunlightTimes}, modified
#'   from \code{amt::time_of_day}. Finds the category of time of day (day vs. night)
#'   given time stamps.
#' @param move A \code{Move*} object or a vector of \code{POSIXct} times.
#' @param crepuscular Should the time of day include crepuscular categories
#'   (dawn & dusk)? Default is FALSE.
#'
#' @return A vector of strings for the time of day.
#' @export
#'
#' @importFrom lubridate %within%
#'
#' @examples
#' NA
#'
#' @import move
time_of_day <- function(move, crepuscular = FALSE) {
  if (!requireNamespace("suncalc", quietly = TRUE)) {
    stop("Please install package `suncalc` first.")
  }
  t <- move@timestamps
  pts <- data.frame(move@coords, as.Date(t))
  names(pts) <- c("lon", "lat", "date")
  sun <- suncalc::getSunlightTimes(data = pts, tz = lubridate::tz(pts$date))

  int.day <- lubridate::interval(sun$sunrise, sun$sunset)
  int.dawn <- lubridate::interval(sun$dawn, sun$sunrise)
  int.dusk <- lubridate::interval(sun$sunset, sun$dusk)

  tod <- c("night", "day")[(t %within% int.day) + 1]

  if (crepuscular) {
    tod[t %within% int.dawn] <- "dawn"
    tod[t %within% int.dusk] <- "dusk"
  }

  if (crepuscular) {
    factor(tod, levels = c("day", "dusk", "night", "dawn"))
  } else {
    factor(tod, levels = c("day", "night"))
  }
}

