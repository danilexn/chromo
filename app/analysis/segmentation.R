# segclust2d alternative functions
tr.segmentation_internal <- function (x, seg.var = NULL, diag.var = NULL, order.var = NULL,
  scale.variable = NULL, Kmax, lmin = NULL, dat = NULL, type = NULL,
  sameSigma = F, subsample_over = 10000, subsample = TRUE,
  subsample_by = NA, ...)
{
  x_nrow <- nrow(x)
  if (missing(Kmax)) {
    Kmax = floor(0.75 * dim(dat)[2]/lmin)
  }
  missing_subsample <- missing(subsample)
  missing_subsample_by <- missing(subsample_by)
  missing_subsample_over <- missing(subsample_over)
  if (subsample) {
    x_nrow <- nrow(x)
    tmp <- subsample(x, subsample_over, subsample_by)
    x <- tmp$x
    subsample_by <- tmp$by
    dat <- dat[, !is.na(x$subsample_ind)]
  }
  else {
    subsample_by <- 1
    x$subsample_ind <- 1:nrow(x)
  }
  if (missing(scale.variable)) {
    scale.variable <- TRUE
  }
  if (scale.variable) {
    dat[1, ] <- scale(dat[1, ])
    dat[2, ] <- scale(dat[2, ])
  }
  else {
    dat[1, ] <- scale(dat[1, ], center = TRUE, scale = FALSE)
    dat[2, ] <- scale(dat[2, ], center = TRUE, scale = FALSE)
  }
  lmin <- base::max(floor(lmin/subsample_by), 2)
  Kmax <- base::min(Kmax, floor(x_nrow/(lmin * subsample_by)))
  if (check_repetition(dat, lmin)) {
    stop("There are repetitions of identical values in the time series larger than lmin, cannot estimate variance for such segment. This is potentially caused by interpolation of missing values or rounding of values.")
  }
  CostLoc <- Gmean_simultanee(dat, lmin = lmin, sameVar = sameSigma)
  res.DynProg <- DynProg_algo_cpp(CostLoc, Kmax)
  output_lavielle <- chooseseg_lavielle(res.DynProg$J.est,
    ...)$Kopt
  res.DynProg$t.est[output_lavielle,]
  return(res.DynProg$t.est[output_lavielle,])
}

tr.segmentation.data.frame <- function (x, Kmax, lmin, type = "home-range", seg.var, diag.var = seg.var, 
  order.var = seg.var[1], coord.names = c("x", "y"), ...)
{
  if (!missing(type)) {
    if (type == "home-range") {
      dat <- t(x[, coord.names])
      seg.var = coord.names
      diag.var = coord.names
      order.var = coord.names[1]
    }
    else if (type == "behavior") {
      if (is.null(seg.var))
        stop("seg.var missing for behavioral segmentation")
      if (length(seg.var) == 1) {
        dat <- t(x[, rep(seg.var, 2)])
      }
      else if (length(seg.var) == 2) {
        dat <- t(x[, seg.var])
      }
      else {
        stop("seg.var must contains either one or two column names")
      }
    }
    else {
      stop("type must be either home-range or behavior")
    }
  }
  else {
    if (!missing(coord.names)) {
      if (missing(seg.var)) {
        seg.var <- coord.names
      }
    }
    if (missing(seg.var)) 
      stop("seg.var missing for behavioral segmentation")
    if (length(seg.var) == 1) {
      dat <- t(x[, rep(seg.var, 2)])
    }
    else if (length(seg.var) == 2) {
      dat <- t(x[, seg.var])
    }
    else {
      stop("seg.var must contains either one or two column names")
    }
  }
  segmented <- tr.segmentation_internal(x, seg.var = seg.var,
    diag.var = diag.var, order.var = order.var, Kmax = Kmax,
    lmin = lmin, dat = dat, type = type, ...)
  return(segmented)
}

segment.lavielle <- function(y, x) {
    outputs <- tr.segmentation.data.frame(data.frame(y, x),lmin = 5,Kmax = 20, seg.var = c("y","x"))
    outputs <- outputs[-which(outputs<=0)]
    ret <- rep(1:length(outputs), c(outputs[1], diff(outputs)))
    return(ret)
}