#' Calculate a home range from a \code{Move} object
#'
#' @description This function currently defaults to an autocorrelated kernel
#'   density estimate, under the assumption that your relocation data is very frequent
#'   (1 hour or less, generally). Several other methods are currently supported,
#'   and more will be added.
#'
#' @param data A \code{Move*} object, or a list of those objects.
#' @param method A home range calculation method:
#'   \itemize{
#'   \item \code{"AKDE"}: Autocorrelated kernel density estimate (default), using the
#'   \code{ctmm} package.
#'   \item \code{"dBBMM"}: Dynamic Brownian Bridge Movement Model from the \code{move}
#'   package.
#'   \item \code{"MCP"}: Minimum convex polygon from the \code{adehabitatHR} package.
#'   \item \code{"KDE"}: Kernel density estimate from the \code{adehabitatHR} package.
#' }
#' @param parallel Should the function be run in parallel? Default is TRUE.
#' @param ... Additional parameters to pass to the home range function of your
#'   choice.
#'
#' @return A home range object.
#' @export
#'
#' @examples
#' NA
#'
#' @import move
home_range <- function(data, method = "AKDE", parallel = TRUE, ...) {
  if (is(data, "MoveStack")) {
    data <- move::split(data)
  }
  if (is.list(data)) {
    if (parallel) {
      requireNamespace("doParallel", quietly = TRUE)
      cores <- parallel::detectCores() - 1
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl, cores = cores)
      opts <- list(preschedule = TRUE)
      parallel::clusterSetRNGStream(cl, 123)
      results <- supressWarnings(plyr::llply(data, home_range, method = method,
                             .progress = "text",
                             .parallel = TRUE,
                             .paropts = list(.options.snow = opts)))
      stopCluster(cl)
    } else {
      results <- plyr::llply(data, home_range, method = method,
                             .progress = "text")
    }
    return(results)
  }
  if (!is(data, "Move")) {
    stop("Data must be a Move (or list of Moves) or MoveStack object
         from the 'Move' package.")
  }
  if (method == "AKDE") {
    #message("Using ctmm::akde()...")
    bird <- as.character(data@idData[["local_identifier"]])
    date <- as.character(data@timestamps[1])
    message(paste0("\n", "Using bird [", bird, "] starting on [", date, "]"))
    #message("Re-projecting CRS to aeqd (Azimuthal Equi-distance) centered on the track...")
    data <- sp::spTransform(data, center=TRUE)
    d <- suppressMessages(ctmm::as.telemetry(data))
    #message("Guesstimating parameters...")
    GUESS <- ctmm::ctmm.guess(d, interactive=FALSE)
    message("Finding best model (this may take a while)...")
    FIT <- ctmm::ctmm.fit(d,GUESS)
    message("Creating home range...")
    UD <- ctmm::akde(d, FIT, ...)
    #message("Done.")
    return(UD)
  }
  if (method == "dBBMM") {
    message("Using move::brownian.bridge.dyn()...")
    message("Creating utilization distribution...")
    dBB <- move::brownian.bridge.dyn(data, ...)
    return(dBB)
  }
  if (method == "MCP") {
    d <- move::move2ade(data)
    MCP <- adehabitatHR::mcp(d, unin = "m", unout = "ha", ...)
    return(MCP)
  }
  if (method == "KDE") {
    d <- move::move2ade(data)
    KDE <- adehabitatHR::kernelUD(d, unin = "m", unout = "ha", ...)
    return(KDE)
  }
}

#' Retrieve the home range size in hectares
#'
#' @param data The output from \code{waterfowl::home_range}.
#' @param size The percentage of the home range to calculate. Default is both
#'   95% (home range) and 50% (core home range).
#' @param method The method that produced the data from \code{waterfowl::home_range}:
#'   \code{"AKDE"}, \code{"dBBMM"}, \code{"MCP"}, \code{"KDE"}.
#'
#' @return A data frame with the home range sizes in hectares.
#' @export
#'
#' @examples
#' NA
#'
#' @import move
home_range_area <- function(data, size = c(95, 50), method = "AKDE") {
  if (is.list(data) && !is(data, "UD")) {
    results <- plyr::llply(data, home_range_area, size = size, method = method,
                           .progress = "text")
    return(results)
  }
  if (method == "AKDE") {
    if (!is(data, "UD")) {
      stop("Data must be a UD object or list of UD objects from the ctmm package.")
    }
    name <- vector(mode = "character", length = length(data))
    df <- data.frame(
      name = data@info$identity
    )
    for (i in 1:length(size)) {
      test <- size[i]
      level <- test / 100
      title <- paste0("akde", test,".ha")
      df[, 1 + i] = (summary(data, level=0.95, level.UD = level,
                             units = FALSE)$CI[[2]]) / 1e4
      colnames(df)[1 + i] <- title
    }
  }
  if (method == "dBBMM") {
    name <- levels(data@DBMvar@dataUnUsedRecords[["local_identifier"]])[1]
    d <- move::getVolumeUD(data)
    dBB50 <- (raster::cellStats(d <= .50, stat = "sum") * raster::res(data)[1]) / 1e4
    dBB95 <- (raster::cellStats(d <= .95, stat = "sum") * raster::res(data)[1]) / 1e4
    df <- data.frame(
      name, dBB50, dBB95
    )
  }
  if (method == "KDE") {
    kde50 <- adehabitatHR::getverticeshr(data, percent = 50, unin = "m", unout = "ha")
    kde95 <- adehabitatHR::getverticeshr(data, percent = 95, unin = "m", unout = "ha" )
    df <- data.frame(
      name = kde50@data[["id"]],
      kde50 = kde50@data[["area"]],
      kde95 = kde95@data[["area"]]
    )
  }
  if (method == "MCP") {
    df <- data.frame(
      name = data@data[["id"]],
      mcp = data@data[["area"]]
    )
  }
  return(df)
}

#' Write a home range to a shapefile
#'
#' @param data The output of \code{waterfowl::home_range}.
#' @param method The method that produced your home range from
#'   \code{waterfowl::home_range}.
#'
#' @return Write a home range shapefile to your current working direcotry.
#' @export
#'
#' @examples
#' NA
#'
#' @import move
write_home_range <- function(data, method = "AKDE") {
  fileName = deparse(substitute(data))
  if (method == "AKDE") {
    ctmm::writeShapefile(
      data, folder = getwd(),
      file = paste0(fileName,"95"),
      level.UD = 0.95, # home range
      level = 0.95 #confidence interval
    )
    ctmm::writeShapefile(
      data, folder = getwd(),
      file = paste0(fileName,"50"),
      level.UD = 0.50, #core home range
      level = 0.95 #confidence interval
    )
  }
  if (method == "dBBM"){
    message("Warning: Only 95% home range implemented.")
    poly = move::raster2contour(data, levels = 0.95)
    raster::shapefile(poly, "dBBMM95.shp")
  }
  if (method == "MCP"){
    stop("MCP not yet implemented.")
  }
  if (method == "KDE"){
    stop("KDE not yet implemented.")
  }
}
