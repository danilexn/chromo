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

calculate.velocity <- function(y, x) {
    vt <- sqrt(diff(y) ^ 2 + diff(x) ^ 2)
    if (length(y) > length(vt)) {
        vt <- vt[1:length(y)]
    }
    return(vt)
}

calculate.velocity.ma <- function(y, x, ma.order) {
    vt <- sqrt(diff(y) ^ 2 + diff(x) ^ 2)
    vt.smooth <- forecast::ma(vt, order = ma.order)
    if (length(y) > length(vt.smooth)) {
        vt.smooth <- vt.smooth[1:length(y)]
    }
    return(vt.smooth)
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
calculate.cumdist <- function(x, y) {
  c(0, cumsum(sp::spDists(as.matrix(cbind(x, y)), segments = TRUE)))
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