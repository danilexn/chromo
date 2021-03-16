# Tronos
# Web App for unsupervised analysis of nuclear oscillations
#
# MIT License
#
# Copyright (c) 2021 Daniel Leon-Perinan (danilexn)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

specCols <- c(
    "#5E4FA2",
    "#3288BD",
    "#66C2A5",
    "#ABDDA4",
    "#E6F598",
    "#FEE08B",
    "#FDAE61",
    "#F46D43",
    "#D53E4F",
    "#9E0142"
)

calculate.spectrum <- function(y, frame) {
    y.spec <-
        spectrum(
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

calculate.freqs <- function(y, frame, sampling.time = 1) {
    y.spec <-
        spectrum(
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

normalize.01 <- function(x, na.rm = FALSE) {
    return((x - min(x)) / (max(x) - min(x)))
}

normalize.z <- function(x, na.rm = FALSE) {
    return(scale(x))
}

normalize.min <- function(x, na.rm = FALSE) {
    return((x) / (min(x)))
}

normalize.max <- function(x, na.rm = FALSE) {
    return((x) / (max(x)))
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
  return(msd.t)
}

segment.timeseries.multiv <-
    function(df,
             grouping,
             particle,
             coords,
             time,
             sampling.time,
             updateProgress = NULL) {
        class.diff <- data.frame()
        if (grouping != "none") {
            groups <- df %>% select(!!sym(grouping)) %>% distinct()
        } else {
            groups <- data.frame(x = 1)
        }
        for (g in groups[, 1]) {
            if (particle != "none") {
                particles <-
                    df %>% filter(!!sym(grouping) == g) %>% select(!!sym(particle)) %>% distinct()
            } else {
                particles <- data.frame(x = 1)
            }
            for (p in particles[, 1]) {
                # Parsing time series
                seq.temp <-
                    df %>% filter(!!sym(grouping) == g) %>% filter(!!sym(particle) == p)

                time.zero <- seq.temp[[time]][1]
                # Selecting columns
                seq.temp.filtered <-
                    as.data.frame(scale(seq.temp[, coords]))
                
                tset <-
                    processTimeseries(
                        ts = seq.temp.filtered,
                        na2zero = TRUE,
                        use.fft = FALSE,
                        dc.trafo = "ash",
                        use.snr = TRUE
                    )
                
                cset <- clusterTimeseries(tset, K = 12)
                segments.calc <-
                    segmentClusters(
                        seq = cset,
                        M = 10,
                        E = 2,
                        nui = 3,
                        S = "icor"
                    )
                segments.loc <- data.frame(segments.calc$segments)
                nsegments <- nrow(segments.calc$segments)
                
                curr.seg <- segments.loc
                names(curr.seg)[2] <- "begin"
                curr.seg <- curr.seg %>% select(begin, end)
                curr.seg$group <- factor(g)
                curr.seg$state <- factor(c(1:nsegments))
                curr.seg$frame.b <- curr.seg$begin
                curr.seg$frame.e <- curr.seg$end
                curr.seg$begin <- curr.seg$begin + time.zero
                curr.seg$end <- curr.seg$end + time.zero
                curr.seg$particle <- factor(p)

                features.segs <- data.frame()
                for (s in 1:nsegments) {
                    beg <- curr.seg[s, "frame.b"]
                    end <- curr.seg[s, "frame.e"]
                    
                    seg.temp <- seq.temp[beg:end, ]
                    seg.temp.mean <-
                        seq.temp.filtered[beg:end, ] %>% summarise_all(mean)
                    
                    spec.s = calculate.spectrum(seg.temp[[coords[1]]], seg.temp[[time]])
                    spec.f = calculate.freqs(seg.temp[[coords[1]]], seg.temp[[time]], sampling.time)
                    
                    spec.s <- spec.s[!is.na(spec.s)]
                    spec.f <- spec.f[!is.na(spec.f)]

                    p1 <- seq.temp[[coords[1]]][beg:end]
                    wave.out <- morlet(p1, 1:length(p1), p2 = 8, dj = 0.1, siglvl = 0.999)
                    # TODO: change the frequency scanning!
                    prd <- wave.out$period[apply(Re(wave.out$wave)[,1:50], 1, which.max)]
                    mod <- lm(prd ~ c(1:length(p1)))
                    mod.slope <- summary(mod)[[4]][2]
                    mod.sign <- mod.slope / abs(mod.slope)

                    features.temp <-
                        data.frame(spec.s, spec.f) %>% top_n(1, spec.s)
                    features.temp <- cbind(features.temp, seg.temp.mean,
                                           data.frame(spec.m = mod.slope, spec.ms = mod.sign))
                    features.segs <- rbind(features.segs, features.temp)
                }
                curr.seg <- cbind(curr.seg, features.segs)
                class.diff <- rbind(class.diff, curr.seg)
                if (is.function(updateProgress)) {
                    text <- paste0("Particle: ", p, ", Group: ", g)
                    updateProgress(text)
                }
            }
        }
        return(class.diff)
    }

segment.timeseries.lavielle <-
    function(df,
             grouping,
             particle,
             coords = c("X", "Y"),
             time,
             K = 20,
             L = 5,
             seg.veloc,
             sampling.time,
             updateProgress = NULL) {
        # TODO: change the velocity calculation to use the own dataframe
        # TODO: change the moving average velocity, to perform it independently from df
        class.diff <- data.frame()
        if (grouping != "none") {
            groups <- df %>% select(!!sym(grouping)) %>% distinct()
        } else {
            groups <- data.frame(x = 1)
        }
        for (g in groups[, 1]) {
            if (particle != "none") {
                particles <-
                    df %>% filter(!!sym(grouping) == g) %>% select(!!sym(particle)) %>% distinct()
            } else {
                particles <- data.frame(x = 1)
            }
            for (p in particles[, 1]) {
                seq.temp <-
                    df %>% filter(!!sym(grouping) == g) %>% filter(!!sym(particle) == p)

                time.zero <- seq.temp[[time]][1]
                if (seg.veloc) {
                    seq.temp$ang.speed <- angular_speed(seq.temp, coord.names = coords)
                    seq.temp$cal.speed <-
                        calculate.velocity(seq.temp[[coords[1]]], seq.temp[[coords[2]]])

                    seq.temp <- seq.temp %>%
                        na.omit()
                    seg.coords <- c("ang.speed", "cal.speed")
                } else {
                    seg.coords <- coords
                }

                seg.1 <-
                    suppressMessages(segmentation(
                        seq.temp,
                        lmin = L,
                        Kmax = K,
                        seg.var = seg.coords
                    ))
                nsegments <- seg.1$Kopt.lavielle
                curr.seg <-
                    cbind(seg.1$outputs[[nsegments]]$segments, seg.1$outputs[[nsegments]]$states[, c(3, 4, 5, 6)])
                curr.seg$group <- factor(g)
                curr.seg$state <-
                    factor(seg.1$outputs[[nsegments]]$states[, c(1)])
                curr.seg$frame.b <- curr.seg$begin
                curr.seg$frame.e <- curr.seg$end
                curr.seg$begin <- curr.seg$begin + time.zero
                curr.seg$end <- curr.seg$end + time.zero
                curr.seg$particle <- factor(p)

                freqs.segs <- data.frame()
                for (s in 1:nsegments) {
                    beg <- seg.1$outputs[[nsegments]]$segments$begin[s]
                    end <- seg.1$outputs[[nsegments]]$segments$end[s]

                    spec.s = calculate.spectrum(seq.temp[[coords[1]]][beg:end], seq.temp[[time]][beg:end])
                    spec.f = calculate.freqs(seq.temp[[coords[1]]][beg:end], seq.temp[[time]][beg:end], sampling.time)

                    spec.s <- spec.s[!is.na(spec.s)]
                    spec.f <- spec.f[!is.na(spec.f)]

                    p1 <- seq.temp[[coords[1]]][beg:end]
                    wave.out <- morlet(p1, 1:length(p1), p2 = 8, dj = 0.1, siglvl = 0.999)
                    # TODO: change the frequency scanning!
                    prd <- wave.out$period[apply(Re(wave.out$wave)[,1:50], 1, which.max)]
                    mod <- lm(prd ~ c(1:length(p1)))
                    mod.slope <- summary(mod)[[4]][2]
                    mod.sign <- mod.slope / abs(mod.slope)
                    freqs.segs <-
                        rbind(freqs.segs,
                              cbind(data.frame(spec.s, spec.f) %>% top_n(1, spec.s),
                              data.frame(spec.m = mod.slope, spec.ms = mod.sign)))
                }
                curr.seg <- cbind(curr.seg, freqs.segs)
                class.diff <- rbind(class.diff, curr.seg)
                if (is.function(updateProgress)) {
                    text <- paste0("Particle: ", p, ", Group: ", g)
                    updateProgress(text)
                }
            }
        }
        return(class.diff)
    }

clusterize.segments <-
    function(class.diff, grouping, particle, segments, updateProgress = NULL) {
        if (grouping != "none" && particle != "none") {
            drop.cols <-
                c("subsample_ind",
                  "subsample_ind.y",
                  "state",
                  "begin",
                  "end",
                  "group",
                  "particle")
            
            class.diff.cluster <-
                class.diff %>% select(-one_of(drop.cols)) %>%
                mutate_each(funs(c(scale(.))))
            
            d_clust <-
                Mclust(
                    as.matrix(class.diff.cluster),
                    G = segments[1]:segments[2],
                    modelNames = mclust.options("emModelNames")
                )
            
            class.diff$cluster <- factor(d_clust$classification)
        } else {
            class.diff$cluster <- factor(class.diff$state)
        }
        return(class.diff)
    }

motifs.per.segment <- function(segments.all, seq.d, seq.q, particle_vars, df_vars_motifs, df_sequence_motifs, motif_window, motif_sample_pct, motif_correlation, updateProgress = NULL) {
    motif.list <- list()
    for (i in 1:nsegments) {
      segments.i <- segments.all %>% filter(cluster == i)
      seq <- c()
      query <- c()
      if (particle_vars != "none") {
          particles <- seq.d %>% select(!!sym(particle_vars)) %>% distinct()
      } else {
          particles <- data.frame(x = 1)
      }
      # TODO: implement progress
      for (p in particles[, 1]) {
        seq.temp <- seq.d %>% filter(!!sym(particle_vars) == p)
        seq.temp <- seq.temp[[df_vars_motifs]]

        segments.i.d <- segments.i %>% filter(group == df_sequence_motifs, particle == p)
        if (nrow(segments.i.d) != 1){
          next
        }
        seq.temp <- seq.temp[segments.i.d$frame.b:segments.i.d$frame.e]

        seq <- c(seq, seq.temp)
      }

      if (particle_vars != "none") {
          particles <- seq.q %>% select(!!sym(particle_vars)) %>% distinct()
      } else {
          particles <- data.frame(x = 1)
      }
      for (p in particles[, 1]) {
        seq.temp <- seq.q %>% filter(!!sym(particle_vars) == p)
        seq.temp <- seq.temp[[df_vars_motifs]]

        segments.i.q <- segments.i %>% filter(group == df_sequence_motifs, particle == p)
        if (nrow(segments.i.q) != 1){
          next
        }
        seq.temp <- seq.temp[segments.i.q$frame.b:segments.i.q$frame.e]

        query <- c(query, seq.temp)
      }

      if (length(seq) * 0.5 < motif_window) {
        next
      }

      text <- paste0("Calculating motifs, segment ", i)
      updateProgress(text)

      motif.seg <-
        motif.discovery(
          seq,
          win = motif_window,
          query,
          pct = motif_sample_pct,
          thr = motif_correlation
        )
      motif.list <- c(motif.list, list(motif.seg))
    }
    return(motif.list)
}

spectrum.significance <- function(df_freqs, particle_vars) {
    m1.r <- lm(spec.s ~ spec.f * sym(particle_vars), data = df_freqs)
    m3.r <- lm(spec.s ~ spec.f, data = df_freqs)
    
    return(anova(m1.r, m3.r)[[6]][2])
}

rotate.axes.maximum <- function(df, coords, particles) {
    df <- df %>% group_by(!!sym(particles)) %>%
        mutate(
            !!sym(coords[1]) := homogeneize.rotation(!!sym(coords[2]),!!sym(coords[1]))[, 2],!!sym(coords[2]) := homogeneize.rotation(!!sym(coords[2]),!!sym(coords[1]))[, 1]
        )
}

homogeneize.rotation <- function(x, y) {
    coords <- cbind(X = x, Y = y)
    rad <- seq(0, pi, l = 20)
    best.rotation <- c()
    for (i in rad) {
        coords.rot <- Rotation(coords, i)
        best.rotation <-
            c(best.rotation, abs(max(coords.rot[, 1]) - min(coords.rot[, 1])))
    }
    
    rot <- rad[which.min(best.rotation)]
    coords.rot <- Rotation(coords, rot)
    return(coords.rot)
}

segmentation.summary <- function(segments) {
    drop.cols <- c("subsample_ind", "subsample_ind.y", "state", "group", "particle")
    segments <- segments %>% select(-one_of(drop.cols))
    return(segments %>% group_by(cluster) %>% summarise_all(mean))
}

segmentation.significance <- function(segments, method, p.correct) {
    ngroups <- segments %>% select(group) %>% distinct()
    if (nrow(ngroups) < 2) {
        return(data.frame(x = "Could not quantify significance. Control group needed."))
    }
    
    mod.1 <-
        glm(group ~ cluster * begin * end,
            data = segments,
            family = "binomial")
    mod.2 <-
        glm(group ~ begin * end, data = segments, family = "binomial")
    mod.3 <-
        glm(group ~ cluster, data = segments, family = "binomial")
    mod.4 <- glm(group ~ 1, data = segments, family = "binomial")

    cluster.durations <- anova(mod.1, mod.2, mod.3, mod.4, test = method)
    # sign.table <-
    #     data.frame(cluster.durations[[ncol(cluster.durations)]][nrow(cluster.durations)],
    #                clusters.durations.intercept[[ncol(clusters.durations.intercept]][nrow(clusters.durations.intercept)],
    #                clusters[[ncol(clusters)]][nrow(clusters)],
    #                durations[[ncol(durations)]][nrow(durations)])
    # sign.table.adj <-
    #     p.adjust(sign.table, method = p.correct, n = length(sign.table))

    test.outputs <- cluster.durations
    return(test.outputs)
}

segmentation.plot <- function(segments, grouping = "none") {
    segno <- segments %>% group_by(group) %>% count(cluster) %>% mutate(freq = n / sum(n))

    p1 <-
        ggplot(data = segments, aes(x = begin, y = cluster, col = cluster)) +
        facet_grid(group ~ .) +
        geom_segment(aes(x = begin, xend = end, yend = cluster), alpha = 0.2) +
        geom_point(aes(x = begin), size = 3, alpha = 0.2) +
        geom_point(aes(x = end), size = 3, alpha = 0.2) +
        theme_bw() +
        xlab("Time") +
        ylab("Cluster") +
        theme(legend.position = "none")

    p2 <- ggplot(data = segno, aes(x = cluster, y = freq, fill = cluster)) +
        facet_grid(group ~ .) +
        geom_col(alpha = 0.5) +
        coord_flip() +
        theme_bw() +
        xlab("Cluster") +
        ylab("Count") +
        theme(legend.position = "none")

    subplot(ggplotly(p1), ggplotly(p2), widths = c(0.7, 0.3), titleX=T)
}

morlet.plot <- function(df, plot.var, sig = 0.95) {
    y <- df[[plot.var]]
    wave.out <-
        morlet(y,
               1:length(y),
               p2 = 8,
               dj = 0.1,
               siglvl = sig)
    wavelet.plot(
        wave.out,
        reverse.y = TRUE,
        key.cols = specCols,
        useRaster = TRUE
    )
}

spectrum.global <-
    function(df,
             segments,
             particle_vars,
             group_vars,
             coord_vars,
             time_vars,
             sampling.time,
             updateProgress = NULL) {
        coord.1 <- coord_vars[1]
        df_freqs <- df %>% group_by(!!sym(particle_vars)) %>%
            mutate(spec.s = calculate.spectrum(!!sym(coord.1),!!sym(time_vars)),
                   spec.f = calculate.freqs(!!sym(coord.1),!!sym(time_vars),sampling.time)) %>%
            na.omit() %>% mutate(spec.s = normalize.z(spec.s))

        df_freqs[[group_vars]] <-
            factor(df_freqs[[group_vars]], levels = sort(unique(df_freqs[[group_vars]])))
        return(df_freqs)
    }

spectrogram.global <-
    function(df,
             particle,
             grouping,
             coord_vars,
             time_vars,
             updateProgress = NULL) {
        coord.1 <- coord_vars[1]

        waves.all <- data.frame()
        if (grouping != "none") {
            groups <- df %>% select(!!sym(grouping)) %>% distinct()
        } else {
            groups <- data.frame(x = 1)
        }
        for (g in groups[, 1]) {
            if (particle != "none") {
                particles <-
                    df %>% filter(!!sym(grouping) == g) %>% select(!!sym(particle)) %>% distinct()
            } else {
                particles <- data.frame(x = 1)
            }
            for (p in particles[, 1]) {
                seq.temp <-
                    df %>% filter(!!sym(grouping) == g) %>% filter(!!sym(particle) == p)
                p1 <- seq.temp[[coord.1]]
                wave.out <- morlet(p1, 1:length(p1), p2 = 8, dj = 0.1)
                Power <- wave.out$Power
                # Select the interval for the morlet spectrogram
                Signif <- as.data.frame(Power[,1:50])
                names(Signif) <- wave.out$period[1:50]

                Signif <- data.frame(t(apply(Signif, 1, function(x) head(names(Signif)[order(-x)],10))))
                Signif <- cbind(Signif, data.frame(frame = seq.temp[[time_vars]]))
                Signif <- Signif %>% pivot_longer(!frame, names_to = "Y", values_to = "spec") %>%
                            mutate_each_(funs(as.numeric), "spec")
                prd.df <- cbind(Signif[,c(1,3)], data.frame(label = seq.temp[[grouping]]))
                waves.all <- rbind(waves.all, prd.df)
                if (is.function(updateProgress)) {
                    text <- paste0("Particle: ", p, ", Group: ", g)
                    updateProgress(text)
                }
            }
        }
        return(waves.all)
    }

spectrogram.per.segment <-
    function(df,segments,
             particle,
             grouping,
             coord_vars,
             time_vars,
             updateProgress = NULL) {
        coord.1 <- coord_vars[1]

        waves.all <- data.frame()
        if (grouping != "none") {
            groups <- df %>% select(!!sym(grouping)) %>% distinct()
        } else {
            groups <- data.frame(x = 1)
        }
        for (g in groups[, 1]) {
            if (particle != "none") {
                particles <-
                    df %>% filter(!!sym(grouping) == g) %>% select(!!sym(particle)) %>% distinct()
            } else {
                particles <- data.frame(x = 1)
            }
            for (p in particles[, 1]) {
                seq.temp <-
                    df %>% filter(!!sym(grouping) == g) %>% filter(!!sym(particle) == p)

                seg.temp <-
                    segments %>% filter(group == g) %>% filter(particle == p)

                if (nrow(seg.temp) < 1) {
                    next
                }

                for (s in 1:nrow(seg.temp)) {
                    beg <- seg.temp$frame.b[s]
                    end <- seg.temp$frame.e[s]

                    seq.temp.seg <- seq.temp[beg:end, ]

                    p1 <- seq.temp.seg[[coord.1]]

                    if (length(p1) < 5) {
                        next
                    }

                    wave.out <- morlet(p1, 1:length(p1), p2 = 8, dj = 0.1)
                    Power <- wave.out$Power
                    # Select the interval for the morlet spectrogram
                    Signif <- as.data.frame(Power[,1:50])
                    names(Signif) <- wave.out$period[1:50]

                    Signif <- data.frame(t(apply(Signif, 1, function(x) head(names(Signif)[order(-x)],10))))
                    Signif <- cbind(Signif, data.frame(frame = seq.temp.seg[[time_vars]]))
                    Signif <- Signif %>% pivot_longer(!frame, names_to = "Y", values_to = "spec") %>%
                                mutate_each_(funs(as.numeric), "spec")
                    cluster <- factor(rep(seg.temp$cluster[s], length(p1)))
                    prd.df <- cbind(Signif[,c(1,3)], data.frame(label = seq.temp.seg[[grouping]], cluster = cluster))
                    waves.all <- rbind(waves.all, prd.df)
                }
                if (is.function(updateProgress)) {
                    text <- paste0("Particle: ", p, ", Group: ", g)
                    updateProgress(text)
                }
            }
        }
        return(waves.all)
    }

g_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <-
        which(sapply(tmp$grobs, function(x)
            x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

spectral.plot <-
    function(df_freqs,
             group_vars,
             sampling.time,
             freq = FALSE) {

        if (freq) {
            p1 <-
                ggplot(df_freqs, aes(spec.f, spec.s), colour = !!sym(group_vars)) +
                xlab("Frequency")
        } else {
            p1 <-
                ggplot(df_freqs, aes((1 / spec.f) / 60, spec.s), colour = !!sym(group_vars)) +
                xlab("Period (minutes)")
        }
        p1 <- p1 + stat_smooth(
                method = "loess",
                span = 0.1,
                se = TRUE,
                aes(
                    color = !!sym(group_vars),
                    fill = !!sym(group_vars)
                ),
                alpha = 0.3
            ) +
            theme_bw() +
            ylab("Normalized spectral density") +
            scale_x_continuous(trans = "log10")

        plt <- ggplotly(p1, height = 500, width = 600)

        return(plt)
    }

spectral.plot.segments <-
    function(df_seg_freqs,
             group_vars,
             sampling.time,
             freq = FALSE) {

        if(is.null(df_seg_freqs)) {
            return(NULL)
        }

        if (freq) {
            p1 <-
                ggplot(df_seg_freqs, aes(spec.f, spec.s), colour = !!sym(group_vars)) +
                xlab("Frequency")
        } else {
            p1 <-
                ggplot(df_seg_freqs, aes((1 / spec.f) / sampling.time, spec.s), colour = !!sym(group_vars)) +
                xlab("Period")
        }
        p1 <- p1 +
            stat_smooth(
                method = "loess",
                span = 0.1,
                se = TRUE,
                aes(color = group, fill = group),
                alpha = 0.3
            ) +
            facet_wrap(cluster ~ ., scales = 'free') +
            theme_bw() +
            ylab("Normalized spectral density") +
            scale_x_continuous(trans = "log10") +
            theme(panel.spacing.x = unit(1, "lines"), panel.spacing.y = unit(4, "lines"))

        plt <- ggplotly(p1, height = 600, width = 600)

        return(plt)
    }

heatmap.plot <-
    function(df_specgram,
             freq = FALSE) {

        if (freq) {
            p1 <-
                df_specgram %>% ggplot(aes(x=frame, y=(1/spec)/60) ) +
                xlab("Frequency")
        } else {
            p1 <-
                df_specgram %>% ggplot(aes(x=frame, y=spec)) +
                xlab("Period")
        }
            p1 <- p1 + stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
                scale_fill_gradientn(colors = specCols) +
                facet_grid(label ~ .) +
                ylab("Period") +
                theme_bw() + theme(legend.position="top",
                                    legend.justification="right",
                                    legend.margin=margin(0,0,0,0),
                                    legend.box.margin=margin(-10,-10,-10,-10))

        plt <- ggplotly(p1)

        return(plt)
    }

heatmap.segmented.plot <-
    function(df_specgram_segmented,
             freq = FALSE) {

        if(is.null(df_specgram_segmented)) {
            return(NULL)
        }

        if (freq) {
            p1 <-
                df_specgram_segmented %>% ggplot(aes(x=frame, y=(1/spec)/60) ) +
                ylab("Frequency")
        } else {
            p1 <-
                df_specgram_segmented %>% ggplot(aes(x=frame, y=spec)) +
                ylab("Period (minutes)")
        }
            p1 <- p1 + stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
                facet_grid(cluster ~ label * .) +
                scale_fill_gradientn(colors = specCols) +
                xlab("Time") +
                theme_bw()

        plt <- ggplotly(p1)

        return(plt)
    }

heatmap.veloc.plot <-
    function(df_vel,
             grouping) {

        names(df_vel)[names(df_vel) == grouping] <- 'label.tronos'

        p1 <- df_vel %>% ggplot(aes(x= frame, y= vel.ma)) +
                xlab("Time") +
                stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
                scale_fill_gradientn(colors = specCols) +
                facet_grid(label.tronos ~ .) +
                ylab("Velocity") +
                theme_bw()

        plt <- ggplotly(p1)

        return(plt)
    }

heatmap.veloc.segmented.plot <-
    function(df_vel_segmented,
             freq = FALSE) {

        if(is.null(df_vel_segmented)) {
            return(NULL)
        }

        p1 <-
            df_vel_segmented %>% ggplot(aes(x=frame, y=vel.ma )) +
            ylab("Velocity") +
            stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
            facet_grid(cluster ~ group * .) +
            scale_fill_gradientn(colors = specCols) +
            xlab("Time") +
            theme_bw()

        plt <- ggplotly(p1)

        return(plt)
    }

spectrum.significance <-
    function(df_freqs,
             segments,
             particle_vars,
             group_vars,
             coord_vars,
             time_vars,
             range) {
        df_freqs$tronos_group <- df_freqs[[group_vars]]
        df_freqs <- df_freqs %>% mutate(spec.f = spec.f / max(spec.f)) %>%
            filter(spec.f < range[2], spec.f > range[1])
        
        m1.r <-
            glm(factor(tronos_group) ~ spec.f * spec.s ,
                data = df_freqs,
                family = "binomial")
        m2.r <-
            glm(factor(tronos_group) ~ spec.f + spec.s ,
                data = df_freqs,
                family = "binomial")
        m3.r <-
            glm(factor(tronos_group) ~ 1 ,
                data = df_freqs,
                family = "binomial")
        anova(m1.r, m2.r, m3.r, test = "Chisq")
    }

spectrum.per.segment <-
    function(df,
             segments,
             coord,
             time_vars,
             grouping,
             particle,
             sampling.time,
             updateProgress = NULL) {
        freqs.segs.grouped <- data.frame()
        if (grouping != "none") {
            groups <- df %>% select(!!sym(grouping)) %>% distinct()
        } else {
            groups <- data.frame(x = 1)
        }
        for (g in groups[, 1]) {
            if (particle != "none") {
                particles <-
                    df %>% filter(!!sym(grouping) == g) %>% select(!!sym(particle)) %>% distinct()
            } else {
                particles <- data.frame(x = 1)
            }
            for (p in particles[, 1]) {
                seq.temp <-
                    df %>% filter(!!sym(grouping) == g) %>% filter(!!sym(particle) == p)

                if (nrow(seq.temp) < 10) {
                    next
                }

                seg.temp <-
                    segments %>% filter(group == g) %>% filter(particle == p)

                if (nrow(seg.temp) < 1) {
                    next
                }

                for (s in 1:nrow(seg.temp)) {
                    beg <- seg.temp$frame.b[s]
                    end <- seg.temp$frame.e[s]
                    
                    seq.temp.seg <- seq.temp[beg:end, ] %>% drop_na(x)
                    if (nrow(seq.temp.seg) < 5) {
                        next
                    }
                    spec.s = calculate.spectrum(seq.temp.seg[[coord]], seq.temp.seg[[time_vars]])
                    spec.f = calculate.freqs(seq.temp.seg[[coord]], seq.temp.seg[[time_vars]], sampling.time)
                    spec.s <- normalize.z(spec.s[!is.na(spec.s)])
                    spec.f <- spec.f[!is.na(spec.f)]
                    
                    cluster <-
                        factor(rep(seg.temp$cluster[s], length(spec.s)))
                    group <- factor(rep(g, length(spec.s)))
                    part <- factor(rep(p, length(spec.s)))
                    
                    freqs.segs.grouped <-
                        rbind(
                            freqs.segs.grouped,
                            data.frame(spec.s, spec.f, cluster, group, part)
                        )
                }
                if (is.function(updateProgress)) {
                    text <- paste0("Particle: ", p, ", Group: ", g)
                    updateProgress(text)
                }
            }
        }
        
        return(freqs.segs.grouped)
    }

velocity.per.segment <-
    function(df,
             segments,
             coords,
             grouping,
             particle,
             vel_ma,
             time_vars,
             updateProgress = NULL) {
        vels.segs.grouped <- data.frame()
        if (grouping != "none") {
            groups <- df %>% select(!!sym(grouping)) %>% distinct()
        } else {
            groups <- data.frame(x = 1)
        }
        for (g in groups[, 1]) {
            if (particle != "none") {
                particles <-
                    df %>% filter(!!sym(grouping) == g) %>% select(!!sym(particle)) %>% distinct()
            } else {
                particles <- data.frame(x = 1)
            }
            for (p in particles[, 1]) {
                seq.temp <-
                    df %>% filter(!!sym(grouping) == g) %>% filter(!!sym(particle) == p)
                if (nrow(seq.temp) < 10) {
                    next
                }
                seg.temp <-
                    segments %>% filter(group == g) %>% filter(particle == p)

                if (nrow(seg.temp) < 1) {
                    next
                }
                for (s in 1:nrow(seg.temp)) {
                    beg <- seg.temp$frame.b[s]
                    end <- seg.temp$frame.e[s]
                    
                    seq.temp.seg <- seq.temp[beg:end, ] %>% drop_na(x)
                    if (nrow(seq.temp.seg) < vel_ma * 2) {
                        next
                    }
                    vel.ma = calculate.velocity.ma(seq.temp.seg[[coords[1]]], seq.temp.seg[[coords[2]]], vel_ma)
                    vel.ma <- vel.ma[!is.na(vel.ma)]
                    
                    cluster <-
                        factor(rep(seg.temp$cluster[s], length(vel.ma)))
                    group <- factor(rep(g, length(vel.ma)))
                    part <- factor(rep(p, length(vel.ma)))
                    diff.ma.t <- length(seq.temp.seg) - length(vel.ma)
                    frame <- seq.temp.seg[[time_vars]][1:(length(seq.temp.seg) - diff.ma.t)]

                    vels.segs.grouped <-
                        rbind(vels.segs.grouped,
                              data.frame(vel.ma, cluster, group, part, frame))
                }
                if (is.function(updateProgress)) {
                    text <- paste0("Particle: ", p, ", Group: ", g)
                    updateProgress(text)
                }
            }
        }
        
        return(vels.segs.grouped)
    }

trajectory.features <-
    function(df, level.vars, ma.order, time, y, segments) {
        # Dataframe reordering for facet plotting
        df_MA <- data.frame(apply(df[, level.vars], 2, SMA, n = ma.order))
        df_MA <-
            df_MA %>% mutate_all(funs((. - min(., na.rm = T)) / (max(., na.rm = T) - min(., na.rm = T))))
        df_MA_long <- melt(df_MA)
        Time = rep(df %>% select(time), length(level.vars))
        
        ### Time series trajectory
        p1 <- ggplot(df, aes(x = time, y = y)) +
            geom_line() +
            xlab("") +
            ylab("Position") +
            theme_bw() +
            geom_vline(
                xintercept = c(best.case$end, best.case$begin),
                color = "black",
                size = 0.5
            )
        
        ### Moving average values of each morphological descriptor and velocity
        p2 <-
            ggplot(data = df_MA_long, aes(
                x = Time,
                y = 1,
                fill = value,
                color = variable
            ))   +
            facet_wrap( ~ variable, ncol = 1, scales = 'free')  +
            scale_fill_gradientn(
                colours = c("#15BFC3", "#FC367A", "white", "gray"),
                na.value = "transparent",
                breaks = c(0, 0.5, 1),
                labels = c(0, 0.5, 1),
                limits = c(0, 1)
            ) +
            theme_bw() +
            ylab('') +
            geom_raster() +
            labs(fill = "Relative\nvalue") +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
            geom_vline(
                xintercept = c(best.case$end, best.case$begin),
                color = "black",
                size = 0.5
            )
        
        p.level <- ggarrange(p1, p2, ncol = 1, heights = c(2, 4))
        return(p.level)
    }

motif.discovery <-
    function(df,
             win = 12,
             q = NULL,
             pct = 1,
             thr = 0.98) {
        motifs <-
            tsmp::analyze(
                df,
                windows = win,
                query = q,
                sample_pct = pct,
                threshold = thr
            )
        return(motifs)
    }

calc.velocities <- function(df, particle, coords, vel_ma) {
    df_vt.smooth <- df %>% group_by(!!sym(particle)) %>%
        mutate(vel.ma = calculate.velocity.ma(!!sym(coords[1]),!!sym(coords[2]), vel_ma)) %>%
        na.omit()
    return(df_vt.smooth)
}

calc.msd <- function(df, particle, grouping, coords, k, updateProgress = NULL) {
    df_msd <- data.frame()
    if (grouping != "none") {
            groups <- df %>% select(!!sym(grouping)) %>% distinct()
    } else {
        groups <- data.frame(x = 1)
    }
    for (g in groups[, 1]) {
        if (particle != "none") {
            particles <-
                df %>% filter(!!sym(grouping) == g) %>% select(!!sym(particle)) %>% distinct()
        } else {
            particles <- data.frame(x = 1)
        }
        for (p in particles[, 1]) {
            seq.temp <-
                    df %>% filter(!!sym(grouping) == g) %>% filter(!!sym(particle) == p)

            seq.temp.x <- seq.temp[[coords[1]]]
            seq.temp.y <- seq.temp[[coords[2]]]

            msd <- calculate.msd(seq.temp.x,seq.temp.y, until = k)
            msd.k <- c(1:k)
            msd.df.temp <- data.frame(msd, msd.k, group = g, particle = p)
            df_msd <- rbind(df_msd, msd.df.temp)
            if (is.function(updateProgress)) {
                text <- paste0("Particle: ", p, ", Group: ", g)
                updateProgress(text)
            }
        }
    }
    return(df_msd)
}

calc.msd.segment <- function(df, segments, particle, grouping, coords, k, updateProgress = NULL) {
    df_msd <- data.frame()
    if (grouping != "none") {
            groups <- df %>% select(!!sym(grouping)) %>% distinct()
    } else {
        groups <- data.frame(x = 1)
    }
    for (g in groups[, 1]) {
        if (particle != "none") {
            particles <-
                df %>% filter(!!sym(grouping) == g) %>% select(!!sym(particle)) %>% distinct()
        } else {
            particles <- data.frame(x = 1)
        }
        for (p in particles[, 1]) {
            seq.temp <-
                    df %>% filter(!!sym(grouping) == g) %>% filter(!!sym(particle) == p)

            if (nrow(seq.temp) < 10) {
                next
            }

            seg.temp <-
                segments %>% filter(group == g) %>% filter(particle == p)

            if (nrow(seg.temp) < 1) {
                next
            }

            for (s in 1:nrow(seg.temp)) {
                beg <- seg.temp$frame.b[s]
                end <- seg.temp$frame.e[s]

                seq.temp.seg <- seq.temp[beg:end, ] %>% drop_na(x)
                if (nrow(seq.temp.seg) < 5) {
                    next
                }

                seq.temp.x <- seq.temp[[coords[1]]]
                seq.temp.y <- seq.temp[[coords[2]]]

                msd <- calculate.msd(seq.temp.x,seq.temp.y, until = k)
                msd.k <- c(1:k)

                cluster <- seg.temp$cluster[s]

                df_msd <-
                    rbind(
                        df_msd,
                        data.frame(msd, msd.k, group = g, particle = p, cluster = factor(cluster))
                    )
            }
            if (is.function(updateProgress)) {
                text <- paste0("Particle: ", p, ", Group: ", g)
                updateProgress(text)
            }
        }
    }
    return(df_msd)
}

significant.velocities <- function(df_vt.smooth, group) {
    wilcox.test(df_vt.smooth$vel.ma ~ df_vt.smooth[[group]], conf.int = TRUE)
}

plot.velocities <-
    function(df_vt.smooth,
             group) {
        p1 <-
            ggplot(df_vt.smooth, aes(
                y = vel.ma,
                fill = !!sym(group),
                x = !!sym(group)
            )) +
            geom_violin(trim = TRUE, alpha = 0.5) +
            geom_boxplot(width = 0.1,
                         trim = TRUE,
                         outlier.shape = NA) +
            theme_bw() +
            ylab("Velocity")

        plt <- ggplotly(p1, autosize = TRUE, height = 500, width = 600)
        # plt <- ggarrange(p1, p2, nrow = 2)
        return(plt)
    }

plot.velocities.segment <-
    function(df_seg_vels) {

        df_seg_vels$group <-
            factor(df_seg_vels$group, levels = sort(levels(df_seg_vels$group), FALSE))

        p2 <-
            ggplot(df_seg_vels, aes(y = vel.ma, fill = group, x = group)) +
            geom_violin(trim = TRUE, alpha = 0.5) +
            geom_boxplot(width = 0.1,
                         trim = TRUE,
                         outlier.shape = NA) +
            facet_wrap(cluster ~ ., ncol = 2) +
            theme_bw() +
            ylab("Velocity")

        plt <- ggplotly(p2, height = 600, width = 600)
        # plt <- ggarrange(p1, p2, nrow = 2)
        return(plt)
    }

plot.msd <-
    function(df_msd) {
        p1 <-
            ggplot(df_msd, aes(
                y = msd,
                color = group,
                fill = group,
                x = msd.k
            )) +
            stat_summary(geom = "line", fun.y = mean) +
            stat_summary(geom = "ribbon", fun.data = mean_cl_normal, alpha = 0.3) +
            theme_bw() +
            ylab("Mean Squared Displacement") +
            xlab("Lag")

        plt <- ggplotly(p1, height = 500, width = 600)
        return(plt)
    }

plot.msd.segment <-
    function(df_msd_segmented) {
        p1 <-
            ggplot(df_msd_segmented, aes(
                y = msd,
                color = group,
                fill = group,
                x = msd.k
            )) +
            stat_summary(geom = "line", fun.y = mean) +
            stat_summary(geom = "ribbon", fun.data = mean_cl_normal, alpha = 0.3) +
            facet_wrap(cluster ~ .) +
            theme_bw() +
            ylab("Mean Squared Displacement") +
            xlab("Lag")

        plt <- ggplotly(p1, height = 500, width = 600)
        return(plt)
    }

clustering.variables <- function(df, variables) {
    # Clustering to retrieve morphology classes
    classes <- Mclust(df %>% select(variables))
    return(classes)
}

clustering.plot <- function(classes) {
    return(plot.Mclust(classes))
}