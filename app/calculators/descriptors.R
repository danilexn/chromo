calculate.spectrum <- function(y, time = NULL) {
    y.spec <-
        stats::spectrum(
            y,
            log = "yes",
            span = 2,
            plot = FALSE,
            na.action = na.remove
        )
    spy <- 2 * y.spec$spec
    if (length(y) > length(spy)) {
        spy <- spy[1:length(y)]
    }
    return(spy)
}

calculate.freqs <- function(y, time = NULL, sampling.time = 1) {
    y.spec <-
        stats::spectrum(
            y,
            log = "yes",
            span = 2,
            plot = FALSE,
            na.action = na.remove
        )
    spx <- y.spec$freq / sampling.time
    if (length(y) > length(spx)) {
        spx <- spx[1:length(y)]
    }
    return(spx)
}

calculate.freqs.max <- function(spx, spy) {
    return(spx[spy == max(spy, na.rm = TRUE)][1])
}

calculate.velocity <- function(df, coords) {
    df_vel <- df %>% select(coords)
    vt <- sqrt(rowSums(as.data.frame(lapply(df_vel, diff, lag=1))^2))
    if (NROW(df_vel) > length(vt)) {
        vt <- vt[1:NROW(df_vel)]
    }
    return(vt)
}

calculate.velocity.ma <- function(df, coords, ma.order) {
    vt <- calculate.velocity(df, coords)
    vt.smooth <- forecast::ma(vt, order = ma.order)
    if (length(vt) > length(vt.smooth)) {
        vt.smooth <- vt.smooth[1:length(vt)]
    }
    return(vt.smooth)
}

calc.angular.speed <- function (df, coords) {
  xx <- diff(df[[coords[1]]])
  if (length(coords) == 1) {
      yy <- diff(df[[coords[1]]])
  } else {
      yy <- diff(df[[coords[2]]])
  }
  b <- sign(xx)
  b[b == 0] <- 1
  bearings <- b * (yy < 0) * pi + atan(xx/yy)
  c(NA, diff(bearings), NA)
}

calculate.ma <- function(x, ma.order) {
    x.smooth <- forecast::ma(x, order = ma.order)
    if (length(x) > length(x.smooth)) {
        x.smooth <- x.smooth[1:length(x)]
    }
    return(x.smooth)
}

# TODO: include the reference
calculate.msd <- function(sx,sy,until=4)
{
  msd.t <- rep(0,until)
  for (dt in 1:until)
  {
    displacement.x <- as.vector(na.omit(sx[(1+dt):length(sx)]) - sx[1:(length(sx)-dt)])
    displacement.y <- as.vector(na.omit(sy[(1+dt):length(sy)]) - sy[1:(length(sy)-dt)])
    sqrdispl <- (displacement.x^2 + displacement.y^2)
    msd.t[dt] <- mean(sqrdispl)
  }
  return(data.frame(msd = msd.t, msd.k = c(1:until)))
}

# TODO: include the reference to amt
calculate.cumdist <- function(df, coords) {
  df_cd <- df %>% select(coords)
  return(c(0, cumsum(sp::spDists(as.matrix(df_cd), segments = TRUE))))
}

spectrum.slope <- function(y) {
    wave.out <- morlet(y, 1:length(y), p2 = 8, dj = 0.1, siglvl = 0.999)
    # TODO: change the frequency scanning!
    prd <- wave.out$period[apply(Re(wave.out$wave)[,1:50], 1, which.max)]
    mod <- lm(prd ~ c(1:length(y)))
    mod.slope <- summary(mod)[[4]][2]
    return(c(mod.slope, mod.slope / abs(mod.slope)))
}

spectrum.global <-
    function(df,
             coords,
             sampling.time,
             updateProgress = NULL) {

        df_freqs <- df %>% group_by(particle) %>%
            mutate(spec.s = calculate.spectrum(!!sym(coords[1])),
                   spec.f = calculate.freqs(!!sym(coords[1]),sampling.time = sampling.time)) %>%
            na.omit() %>% mutate(spec.s = normalize.z(spec.s))

        return(as.data.frame(df_freqs))
    }

kurtosis <- function (x, na.rm = FALSE)
{
    if (is.matrix(x))
        apply(x, 2, kurtosis, na.rm = na.rm)
    else if (is.vector(x)) {
        if (na.rm)
            x <- x[!is.na(x)]
        n <- length(x)
        n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
    }
    else if (is.data.frame(x))
        sapply(x, kurtosis, na.rm = na.rm)
    else kurtosis(as.vector(x), na.rm = na.rm)
}

skewness <- function (x, na.rm = FALSE)
{
    if (is.matrix(x))
        apply(x, 2, skewness, na.rm = na.rm)
    else if (is.vector(x)) {
        if (na.rm)
            x <- x[!is.na(x)]
        n <- length(x)
        (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
    }
    else if (is.data.frame(x))
        sapply(x, skewness, na.rm = na.rm)
    else skewness(as.vector(x), na.rm = na.rm)
}