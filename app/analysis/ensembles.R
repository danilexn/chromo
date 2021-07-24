calc.velocities <- function(df, coords, vel_ma) {
    df_vt.smooth <- df %>% group_by(group, particle) %>%
        mutate(vel.ma = calculate.velocity.ma(., coords, vel_ma)) %>%
        na.omit()
    return(as.data.frame(df_vt.smooth))
}

calc.densities <- function(df, coords, vel_ma) {
    df_vt.smooth <- df %>% group_by(group, particle) %>%
        mutate(vel.ma = calculate.ma(!!sym(coords[1]), vel_ma)) %>%
        na.omit()
    return(as.data.frame(df_vt.smooth))
}

calc.msd <- function(df, coords, k) {
    df_msd <- df %>% group_by(group, particle) %>%
        mutate(msd.x = !!sym(coords[1]), msd.y = !!sym(coords[2])) %>%
        do(possibly(calculate.msd, otherwise = data.frame(msd = NA, msd.k = NA))(.$msd.x,.$msd.y, k)) %>%
        na.omit()
    return(df_msd)
}

calc.msd.segment <- function(df, coords, k) {
    df_msd <- df %>% group_by(group, particle, cluster) %>%
        mutate(msd.x = !!sym(coords[1]), msd.y = !!sym(coords[2])) %>%
        do(possibly(calculate.msd, otherwise = data.frame(msd = NA, msd.k = NA))(.$msd.x,.$msd.y, k)) %>%
        na.omit()
    return(df_msd)
}

calc.cor <- function(df, vars, sel_group) {
    df_group <- df[df$group == sel_group, vars]
    return(df_group)
}

spectrogram.calculate <- function(y, frame) {
    wave.out <- morlet(y, 1:length(y), p2 = 8, dj = 0.1)
    Power <- wave.out$Power
    # TODO: Select the interval for the morlet spectrogram
    Signif <- as.data.frame(Power[,1:50])
    names(Signif) <- wave.out$period[1:50]

    # TODO: speedup this with Rcpp, maybe (?)
    Signif <- data.frame(t(apply(Signif, 1, function(x) head(names(Signif)[order(-x)],10))))
    Signif <- cbind(Signif, data.frame(frame = frame))
    Signif <- Signif %>% pivot_longer(!frame, names_to = "Y", values_to = "spec") %>%
                mutate_each_(funs(as.numeric), "spec")
    return(as.data.frame(Signif))
}
