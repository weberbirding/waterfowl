% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/landcover.R
\name{landcover_points}
\alias{landcover_points}
\title{Determine land cover class at specified point locations}
\usage{
landcover_points(points, raster.source = "FedData", add = FALSE)
}
\arguments{
\item{points}{Spatial data of class \code{Move*}, \code{SpatialPointsDataFrame},
or \code{SpatVector}, or a list of those classes.}

\item{raster.source}{A raster. By default, downloads the National Land Cover
Database using the \code{FedData} package. Otherwise, accepts any raster that can
work with \code{terra::rast()}.}

\item{add}{Should the output be added to the input spatial data? If FALSE (default),
returns a one-column data frame. If TRUE, returns the input object with an
additional data column "landcover".}
}
\value{
A data frame of the land cover types, or the original data with an
added column of land cover types if add = TRUE.
}
\description{
Determine land cover class at specified point locations
}
\examples{
NA

}
