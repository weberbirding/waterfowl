#' Calculate land cover class percentages for a polygon
#'
#' @param data Spatial data of class \code{Move*} or \code{SpatialPointsDataFrame},
#'   or a list of those classes.
#' @param raster.source A raster. By default, downloads the National Land Cover
#'   Database using the \code{FedData} package. Otherwise, accepts any raster that can
#'   work with \code{terra::rast()}.
#' @param buffer An optional buffer to add to your spatial data, in meters.
#' @param collapse If TRUE and the input was a list, collapses the list into a
#'   single data frame.
#' @param pivot_wide If TRUE, pivots the resulting data frame so that each value
#'   type (by default land cover classes) is in its own column. Useful if you intend
#'   to model habitat selection.
#'
#' @return A data frame of percentages.
#' @export
#'
#' @examples
#' NA
#' @importFrom magrittr %>%
#' @importFrom methods is
landcover_percent <- function(data, raster.source = "FedData",
                          buffer = 0, collapse = FALSE, pivot_wide = TRUE){
  value <- coverage_area <- . <- nlcd_legend <- Class <- percent_cover <-
    Class <- NULL
  if (is.list(data)){
    results <- plyr::llply(data, landcover_percent, raster.source = raster.source,
                           buffer = buffer, collapse = collapse,
                           pivot_wide = pivot_wide, .progress = "text")
    return(results)
  }
  if (is(raster.source, "character")){
    raster <- FedData::get_nlcd(
      template = data,
      label = "landcover",
      dataset = "landcover",
      force.redo = TRUE
    )
  } else{
    raster <- raster.source
  }
  if (!is(raster, "SpatRaster")){
    spat <- terra::rast(raster)
  } else {
    spat <- raster
  }
  if (is(data, "Move") || is(data, "MoveStack")){
    ade <- move::move2ade(data)
    point_names <- rownames(ade@coords)
    data <- sp::spTransform(ade, terra::crs(spat))
  } else if (is(data, "SpatialPointsDataFrame")){
    point_names <- rownames(data@coords)
    data <- sp::spTransform(data, terra::crs(spat))
  } else {
    stop("Data must be of class Move, MoveStack or SpatialPointsDataFrame")
  }
  if (buffer > 0){
    sf_data <- sf::st_as_sf(data, coords = c("location_long", "location_lat"),
                        crs = terra::crs(data), agr = "identity")
    sf_data <- sf::st_transform(sf_data, terra::crs(raster))
    rownames(sf_data) <- point_names
    data <- sf::st_buffer(sf_data, dist = buffer)
  }
  message("Extracting raster values...")
  area <- exactextractr::exact_extract(spat, data, coverage_area = TRUE, progress = TRUE)
  get_area <- function(area){
    area.total <- area %>%
      dplyr::group_by(value) %>%
      dplyr::summarise(
        total_area_ha = sum(coverage_area)/10000,
        percent_cover = sum(coverage_area)/sum(area$coverage_area)*100
      ) %>%
      dplyr::rename(ID = value) %>%
      dplyr::left_join(., nlcd_legend, by = "ID") %>%
      dplyr::select(
        Class,
        #total_area_ha,
        percent_cover
      )
    return(area.total)
  }
  if (is.list(area)){
    message("Summarizing extracted values...")
    area.list <- plyr::llply(area, get_area, .progress = "text")
    names(area.list) <- point_names
    if (collapse == TRUE){
      area <- dplyr::bind_rows(area.list, .id = "point_id")
    }
  } else {
    area <- get_area(area)
  }
  if (pivot_wide == TRUE){
    area <- area %>%
      tidyr::pivot_wider(names_from = Class, values_from = percent_cover, values_fill = 0L)
  }
  return(as.data.frame(area))
}

#' Determine land cover class at specified point locations
#'
#' @param points Spatial data of class \code{Move*}, \code{SpatialPointsDataFrame},
#'   or \code{SpatVector}, or a list of those classes.
#' @param raster.source A raster. By default, downloads the National Land Cover
#'   Database using the \code{FedData} package. Otherwise, accepts any raster that can
#'   work with \code{terra::rast()}.
#' @param df Should the output be a data frame? If FALSE, returns the input object
#'   with an additional data column "landcover".
#'
#' @return The original data with an added column of land cover types, or a data
#'   frame of the land cover types if df = TRUE.
#' @export
#'
#' @examples
#' NA
#'
#' @importFrom magrittr %>%
#' @importFrom methods is
landcover_points <- function(points, raster.source = "FedData", df = FALSE){
  tag_id <- local_identifier <- location_long <- location_lat <- timestamp <-
    NULL
  if (is.list(points)){
    results <- plyr::llply(points, landcover_points, raster.source = raster.source,
                           df = df, .progress = "text")
    return(results)
  }
  if (is(raster.source, "character")){
    raster <- FedData::get_nlcd(
      template = vec,
      label = "landcover",
      dataset = "landcover",
      force.redo = TRUE
    )
  } else{
    raster <- raster.source
  }
  if (is(raster, "SpatRaster")){
    spat <- terra::rast(raster)
  } else {
    spat <- raster
  }
  if (is(points, "Move") || is(points, "MoveStack")){
    ade <- move::move2ade(points) %>%
      sp::spTransform(ade, terra::crs(spat))
    vec <- terra::vect(ade)
    names <- rownames(ade@coords)
  } else if (is(points, "SpatialPointsDataFrame")){
    ade <- sp::spTransform(points, terra::crs(spat))
    vec <- terra::vect(ade)
    names <- rownames(points@coords)
  } else if (is(points, "SpatVector")){
    vec <- terra::project(points, terra::crs(spat))
  } else {
    stop("Data must be of class Move, SpatialPointsDataFrame, or SpatVector")
  }
  landcover <- terra::extract(spat, vec, bind = TRUE)
  if (df == FALSE){
    if (is(points, "SpatVector")){
      return(landcover)
    } else {
      points@data$landcover <- as.character(landcover[[2]][,1])
      return(points)
    }
  } else if (df == TRUE){
    if (is(points, "Move") || is(points, "MoveStack")){
      points@data$landcover <- as.character(landcover[[2]][,1])
      df <- as.data.frame(points) %>%
        dplyr::select(
          tag_id,
          local_identifier,
          landcover,
          location_long,
          location_lat,
          timestamp
        )
    } else if (is(points, "SpatialPointsDataFrame")){
      points@data$landcover <- as.character(landcover[[2]][,1])
      df <- as.data.frame(points)
    } else if (is(points, "SpatVector")){
      vec.df <- terra::project(landcover, terra::crs(points))
      df <- as.data.frame(vec.df)
    }
    return(df)
  }
}
