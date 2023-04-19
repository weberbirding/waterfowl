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
#'
#' @import move
#' @importFrom methods is
landcover_percent <- function(data, raster.source = "FedData",
                          buffer = 0, collapse = FALSE, pivot_wide = TRUE){
  if (is.list(data)){
    results <- plyr::llply(data, landcover_percent, raster.source = raster.source,
                           buffer = buffer, collapse = collapse,
                           pivot_wide = pivot_wide, .progress = "text")
    return(results)
  }
  if (is.na(terra::crs(points))) {
    warning("CRS is NA. Assuming WGS84 unprojected lon/lat.")
    terra::crs(points) <- "+proj=longlat +datum=WGS84 +no_defs"
  }
  if (is(raster.source, "character")){
    # create template with slight buffer to avoid extracting from a boundary
    box <- sf::st_bbox(points)
    extent <- box + c(-1, -1, 1, 1)
    extent_points <- data.frame(long = c(extent[1], extent[3]),
                                lat = c(extent[2], extent[4]))
    spdf <- sp::SpatialPoints(extent_points, proj4string = terra::crs(points))
    # grab landcover raster clipped to extent
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
  nlcd_legend <- FedData::nlcd_colors()
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
#' @param add Should the output be added to the input spatial data? If FALSE (default),
#'   returns a one-column data frame. If TRUE, returns the input object with an
#'   additional data column "landcover".
#'
#' @return A data frame of the land cover types, or the original data with an
#'   added column of land cover types if add = TRUE.
#' @export
#'
#' @examples
#' NA
#'
#' @import move
#' @importFrom methods is
landcover_percent <- function(data, raster.source = "FedData",
                             buffer = 0, pivot_wide = TRUE){
  if (is.list(data)){
    results <- plyr::llply(data, landcover_percent, raster.source = raster.source,
                           buffer = buffer, collapse = collapse,
                           pivot_wide = pivot_wide, .progress = "text")
    return(results)
  }
  if (is.na(terra::crs(data))) {
    warning("CRS is NA. Assuming WGS84 unprojected lon/lat.")
    terra::crs(data) <- "+proj=longlat +datum=WGS84 +no_defs"
  }
  if (is(raster.source, "character")){
    # create template with slight buffer to avoid extracting from a boundary
    box <- sf::st_bbox(data)
    extent <- box + c(-1, -1, 1, 1)
    extent_points <- data.frame(long = c(extent[1], extent[3]),
                                lat = c(extent[2], extent[4]))
    spdf <- sp::SpatialPoints(extent_points, proj4string = terra::crs(data))
    # grab landcover raster clipped to extent
    raster <- FedData::get_nlcd(
      template = spdf,
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
  nlcd_legend <- FedData::nlcd_colors()
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
    area <- dplyr::bind_rows(area.list, .id = "point_id")
  } else {
    area <- get_area(area)
  }
  if (pivot_wide == TRUE){
    area <- area %>%
      tidyr::pivot_wider(names_from = Class,
                         values_from = percent_cover,
                         values_fill = 0L)
  }
  return(as.data.frame(area))
}

## quiets concerns of R CMD check for variables that appear in pipelines
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(".", "value", "coverage_area",
                           "nlcd_legend", "Class",
                           "percent_cover", "tag_id",
                           "local_identifier", "location_long",
                           "location_lat", "timestamp"))
}
