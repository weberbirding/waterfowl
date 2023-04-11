#' Pivot Ornitela accelerometry data to long format (xyzxyzxyz...)
#'
#' @description This function takes a data frame made from a .csv downloaded
#'   directly from the Ornitela website, and pivots it into the format
#'   "xyzxyzxyz...", which is used for analyses in the package \code{rabc} on
#'   Github.
#'
#' @param data A data frame of an Ornitela .csv for acc data ("SENSORS").
#' @param burstMinApart The frequency of ACC bursts in minutes.
#' @param pointsPerBurst The number of points in a burst, a.k.a. hertz * seconds.
#' @param removeIncompleteBursts Should incomplete bursts (bursts with fewer
#'   than \code{pointsPerBurst}) be removed from the final data frame? The default
#'   is \code{TRUE}.
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
    dplyr::filter(., .$datatype == "SENSORS") %>%
    dplyr::select(device_id, UTC_datetime, acc_x, acc_y, acc_z)
  data$UTC_datetime <- as.POSIXct(data$UTC_datetime,
                                  tryFormats = c("%m/%d/%y %H:%M:%S",
                                                 "%m/%d/%y %H:%M",
                                                 "%m/%d/%Y %H:%M",
                                                 "%m/%d/%Y %H:%M:%S"), tz = "GMT")
  uniqueBursts <- unique(data$UTC_datetime)
  bursts <- c(uniqueBursts[1])
  for (i in 2:length(uniqueBursts)){
    if (base::difftime(uniqueBursts[i],
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
    result <- result[stats::complete.cases(result),]
  }
  return(result)
}

calculate_feature_extra <- function(df_raw = NULL, winlen_dba = 9, axis_num = 3)
{
  if (is.null(df_raw)) {
    stop("Please provide a valid data.frame!")
  }
  if (!exists("winlen_dba")) {
    stop("Provide a valid window length for running mean")
  }
  if (!is.numeric(winlen_dba)) {
    stop("Winlen_dba should be a number (less than your burst hertz)")
  }
  if (axis_num == 3 & (ncol(df_raw) - 1)%%3 != 0) {
    stop("Number of data column (not including the label column) should be three multiples.")
  }
  if (axis_num == 2 & (ncol(df_raw) - 1)%%2 != 0) {
    stop("Number of data column (not including the label column) should be two multiples.")
  }
  if (is.unsorted(as.character(df_raw[, ncol(df_raw)]))) {
    warning("This function expects the input data sorted by function order_acc.")
  }
  row_num <- dim(df_raw)[[1]]
  col_num <- dim(df_raw)[[2]]
  val_range <- max(as.numeric(as.matrix(df_raw[, -col_num]))) -
    min(as.numeric(as.matrix(df_raw[, -col_num])))
  df_raw <- as.data.frame(df_raw)
  if (axis_num == 3) {
    sub_x <- df_raw[, seq(from = 1, to = col_num - 1, by = 3)]
    sub_y <- df_raw[, seq(from = 2, to = col_num - 1, by = 3)]
    sub_z <- df_raw[, seq(from = 3, to = col_num - 1, by = 3)]
    matsub_x <- as.matrix(sub_x)
    matsub_y <- as.matrix(sub_y)
    matsub_z <- as.matrix(sub_z)
    matmsa <- abs(sqrt((matsub_x^2 + matsub_y^2 + matsub_z^2)) - 1000)
    msa <- rowMeans(matmsa)
    matpitch <- atan2(matsub_y, sqrt(matsub_x^2 + matsub_z^2)) * (180 / pi)
    pitch <- rowMeans(matpitch)
    matroll <- atan2(matsub_x, sqrt(matsub_y^2 + matsub_z^2)) * (180 / pi)
    roll <- rowMeans(matroll)
    matyaw <- atan2(matsub_z, sqrt(matsub_x^2 + matsub_y^2)) * (180 / pi)
    yaw <- rowMeans(matyaw)
    mat_static_x <- t(caTools::runmean(x = t(matsub_x), k = winlen_dba,
                                       alg = "C", endrule = "mean"))
    mat_static_y <- t(caTools::runmean(x = t(matsub_y), k = winlen_dba,
                                       alg = "C", endrule = "mean"))
    mat_static_z <- t(caTools::runmean(x = t(matsub_z), k = winlen_dba,
                                       alg = "C", endrule = "mean"))
    mat_dynamic_x <- matsub_x - mat_static_x
    mat_dynamic_y <- matsub_y - mat_static_y
    mat_dynamic_z <- matsub_z - mat_static_z
    static_x <- rowSums(abs(mat_static_x))
    static_y <- rowSums(abs(mat_static_y))
    static_z <- rowSums(abs(mat_static_z))
    dynamic_x <- rowSums(abs(mat_dynamic_x))
    dynamic_y <- rowSums(abs(mat_dynamic_y))
    dynamic_z <- rowSums(abs(mat_dynamic_z))
    odba <- rowMeans(abs(mat_dynamic_x) + abs(mat_dynamic_y) + abs(mat_dynamic_z))
    vedba <- rowMeans(sqrt((mat_dynamic_x^2) + (mat_dynamic_y^2) + (mat_dynamic_z^2)))
    mean.dxy <- rowMeans(matsub_x - matsub_y)
    mean.dyz <- rowMeans(matsub_y - matsub_z)
    mean.dxz <- rowMeans(matsub_x - matsub_z)
    sd.dxy <- apply((mat_dynamic_x - mat_dynamic_y), 1, stats::sd)
    sd.dyz <- apply((mat_dynamic_y - mat_dynamic_z), 1, stats::sd)
    sd.dxz <- apply((mat_dynamic_x - mat_dynamic_z), 1, stats::sd)
    sk.dx <- apply(mat_dynamic_x, 1, e1071::skewness)
    ku.dx <- apply(mat_dynamic_x, 1, e1071::kurtosis)
    sk.dy <- apply(mat_dynamic_y, 1, e1071::skewness)
    ku.dy <- apply(mat_dynamic_y, 1, e1071::kurtosis)
    sk.dz <- apply(mat_dynamic_z, 1, e1071::skewness)
    ku.dz <- apply(mat_dynamic_z, 1, e1071::kurtosis)
    cov.dxy <- sapply(1:nrow(mat_dynamic_x), function(i) {stats::cov(mat_dynamic_x[i,], mat_dynamic_y[i,])})
    cov.dxz <- sapply(1:nrow(mat_dynamic_x), function(i) {stats::cov(mat_dynamic_x[i,], mat_dynamic_z[i,])})
    cov.dyz <- sapply(1:nrow(mat_dynamic_y), function(i) {stats::cov(mat_dynamic_y[i,], mat_dynamic_z[i,])})
    cor.dxy <- sapply(1:nrow(mat_dynamic_x), function(i) {stats::cor(mat_dynamic_x[i,], mat_dynamic_y[i,])})
    cor.dxz <- sapply(1:nrow(mat_dynamic_x), function(i) {stats::cor(mat_dynamic_x[i,], mat_dynamic_z[i,])})
    cor.dyz <- sapply(1:nrow(mat_dynamic_y), function(i) {stats::cor(mat_dynamic_y[i,], mat_dynamic_z[i,])})
    quant.dx <- apply(mat_dynamic_x, 1, stats::quantile)
    quant.dy <- apply(mat_dynamic_y, 1, stats::quantile)
    quant.dz <- apply(mat_dynamic_z, 1, stats::quantile)
    max.dx <- apply(mat_dynamic_x, 1, max)
    min.dx <- apply(mat_dynamic_x, 1, min)
    amp.dx <- abs(max.dx-min.dx)
    max.dy <- apply(mat_dynamic_y, 1, max)
    min.dy <- apply(mat_dynamic_y, 1, min)
    amp.dy <- abs(max.dy-min.dy)
    max.dz <- apply(mat_dynamic_z, 1, max)
    min.dz <- apply(mat_dynamic_z, 1, min)
    amp.dz <- abs(max.dz-min.dz)
    # summary statistics from accelerateR package
    mdocp_x <- sapply(1:nrow(mat_dynamic_x),
                      function(i){mean(mat_dynamic_x[i,-1] -
                                         mat_dynamic_x[i,-length(mat_dynamic_x[i,])])})
    mdocp_y <- sapply(1:nrow(mat_dynamic_y),
                      function(i){mean(mat_dynamic_y[i,-1] -
                                         mat_dynamic_y[i,-length(mat_dynamic_y[i,])])})
    mdocp_z <- sapply(1:nrow(mat_dynamic_z),
                      function(i){mean(mat_dynamic_z[i,-1] -
                                         mat_dynamic_z[i,-length(mat_dynamic_z[i,])])})
    sdocp_x <- sapply(1:nrow(mat_dynamic_x),
                      function(i){stats::sd(mat_dynamic_x[i,-1] -
                                       mat_dynamic_x[i,-length(mat_dynamic_x[i,])])})
    sdocp_y <- sapply(1:nrow(mat_dynamic_y),
                      function(i){stats::sd(mat_dynamic_y[i,-1] -
                                       mat_dynamic_y[i,-length(mat_dynamic_y[i,])])})
    sdocp_z <- sapply(1:nrow(mat_dynamic_z),
                      function(i){stats::sd(mat_dynamic_z[i,-1] -
                                       mat_dynamic_z[i,-length(mat_dynamic_z[i,])])})
    x_mean <- rowMeans(matsub_x)
    y_mean <- rowMeans(matsub_y)
    z_mean <- rowMeans(matsub_z)
    x_sd <- apply(matsub_x, 1, stats::sd)
    y_sd <- apply(matsub_y, 1, stats::sd)
    z_sd <- apply(matsub_z, 1, stats::sd)
    x_cv <- x_sd/x_mean
    y_cv <- y_sd/y_mean
    z_cv <- z_sd/z_mean
    x_icv <- x_mean/x_sd
    y_icv <- y_mean/y_sd
    z_icv <- z_mean/z_sd
    xyz_sqrt_ss <- sqrt(x_mean^2 + y_mean^2 + z_mean^2)
    # from AcceleRater
    lines_cross <- function(x, y) {
      cross <- 0
      for (i in 2:length(x)) {
        if(x[i-1] > y[i-1] && x[i] < y[i]) {
          cross <- cross + 1
        } else if(x[i-1] < y[i-1] && x[i] > y[i]) {
          cross <- cross + 1
        }
      }
      return(cross)
    }
    crossings.xy <- sapply(1:nrow(matsub_x),
                           function(i){lines_cross(x = matsub_x[i,], y = matsub_y[i,])})
    crossings.xz <- sapply(1:nrow(matsub_x),
                           function(i){lines_cross(x = matsub_x[i,], y = matsub_z[i,])})
    crossings.yz <- sapply(1:nrow(matsub_y),
                           function(i){lines_cross(x = matsub_y[i,], y = matsub_z[i,])})
    norm.x <- apply(matsub_x, 1, base::norm, type = "2")
    norm.y <- apply(matsub_y, 1, base::norm, type = "2")
    norm.z <- apply(matsub_z, 1, base::norm, type = "2")
    df_feature <- data.frame(msa, pitch, roll, yaw,
                             static_x, static_y, static_z,
                             dynamic_x, dynamic_y, dynamic_z,
                             odba,
                             vedba,
                             mean.dxy,mean.dyz,mean.dxz,
                             sd.dxy,sd.dyz,sd.dxz,
                             cov.dxy,cov.dxz,cov.dyz,
                             cor.dxy,cor.dxz,cor.dyz,
                             sk.dx,ku.dx,
                             sk.dy,ku.dy,
                             sk.dz,ku.dz,
                             pct25.dx=quant.dx[2,],pct50.dx=quant.dx[3,],pct75.dx=quant.dx[4,],
                             pct25.dy=quant.dy[2,],pct50.dy=quant.dy[3,],pct75.dy=quant.dy[4,],
                             pct25.dz=quant.dz[2,],pct50.dz=quant.dz[3,],pct75.dz=quant.dz[4,],
                             max.dx,min.dx, amp.dx,
                             max.dy,min.dy, amp.dy,
                             max.dz,min.dz, amp.dz,
                             mdocp_x, mdocp_y, mdocp_z,
                             sdocp_x, sdocp_y, sdocp_z,
                             x_cv, y_cv, z_cv,
                             x_icv, y_icv, z_icv,
                             xyz_sqrt_ss,
                             crossings.xy, crossings.xz, crossings.yz,
                             norm.x, norm.y, norm.z)
  }
  else if (axis_num == 2) {
    sub_x <- df_raw[, seq(from = 1, to = col_num - 1, by = 3)]
    sub_y <- df_raw[, seq(from = 2, to = col_num - 1, by = 3)]
    matsub_x <- as.matrix(sub_x)
    matsub_y <- as.matrix(sub_y)
    mat_static_x <- t(caTools::runmean(x = t(matsub_x), k = winlen_dba, alg="C", endrule="mean"))
    mat_static_y <- t(caTools::runmean(x = t(matsub_y), k = winlen_dba, alg="C", endrule="mean"))
    mat_dynamic_x <- matsub_x - mat_static_x
    mat_dynamic_y <- matsub_y - mat_static_y
    static_x <-rowSums(abs(mat_static_x))
    static_y <-rowSums(abs(mat_static_y))
    dynamic_x <- rowSums(abs(mat_dynamic_x))
    dynamic_y <- rowSums(abs(mat_dynamic_y))
    mean.dxy <- rowMeans(matsub_x-matsub_y)
    sd.dxy <- apply((mat_dynamic_x-mat_dynamic_y), 1, stats::sd)
    sk.dx <- apply(mat_dynamic_x, 1, e1071::skewness)
    ku.dx <- apply(mat_dynamic_x, 1, e1071::kurtosis)
    sk.dy <- apply(mat_dynamic_y, 1, e1071::skewness)
    ku.dy <- apply(mat_dynamic_y, 1, e1071::kurtosis)
    cov.dxy <- sapply(1:nrow(mat_dynamic_x),function(i){stats::cov(mat_dynamic_x[i,], mat_dynamic_y[i,])})
    cor.dxy <- sapply(1:nrow(mat_dynamic_x),function(i){stats::cor(mat_dynamic_x[i,], mat_dynamic_y[i,])})
    quant.dx <- apply(mat_dynamic_x, 1, stats::quantile)
    quant.dy <- apply(mat_dynamic_y, 1, stats::quantile)
    max.dx <- apply(mat_dynamic_x, 1, max)
    min.dx <- apply(mat_dynamic_x, 1, min)
    amp.dx <- abs(max.dx-min.dx)
    max.dy <- apply(mat_dynamic_y, 1, max)
    min.dy <- apply(mat_dynamic_y, 1, min)
    amp.dy <- abs(max.dy-min.dy)
    df_feature <- data.frame(static_x, static_y, dynamic_x, dynamic_y,
                             mean.dxy, sd.dxy, cov.dxy, cor.dxy, sk.dx,ku.dx, sk.dy,ku.dy,
                             pct25.dx=quant.dx[2,],pct50.dx=quant.dx[3,],pct75.dx=quant.dx[4,],
                             pct25.dy=quant.dy[2,],pct50.dy=quant.dy[3,],pct75.dy=quant.dy[4,],
                             max.dx,min.dx, amp.dx, max.dy,min.dy, amp.dy)
  }
  else if (axis_num == 1) {
    sub_x <- df_raw[, seq(from = 1, to = col_num - 1, by = 3)]
    matsub_x <- as.matrix(sub_x)
    mat_static_x <- t(caTools::runmean(x = t(matsub_x), k = winlen_dba, alg="C", endrule="mean"))
    mat_dynamic_x <- matsub_x - mat_static_x
    static_x <-rowSums(abs(mat_static_x))
    dynamic_x <- rowSums(abs(mat_dynamic_x))
    sk.dx <- apply(mat_dynamic_x, 1, e1071::skewness)
    ku.dx <- apply(mat_dynamic_x, 1, e1071::kurtosis)
    quant.dx <- apply(mat_dynamic_x, 1, stats::quantile)
    max.dx <- apply(mat_dynamic_x, 1, max)
    min.dx <- apply(mat_dynamic_x, 1, min)
    amp.dx <- abs(max.dx-min.dx)
    df_feature <- data.frame(static_x, dynamic_x,sk.dx, ku.dx,
                             pct25.dx=quant.dx[2,],pct50.dx=quant.dx[3,],pct75.dx=quant.dx[4,],
                             max.dx,min.dx, amp.dx)
  }
  else {
    stop("Please provide valid number of axis from 1, 2, or 3.")
  }
  return(df_feature)
}
