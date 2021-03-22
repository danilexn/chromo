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

create.df.clusters <- function(df_cluster, df_segment) {
    df_cluster <- df_cluster %>% arrange(group, particle, segmenttronos)
    df_segment <- df_segment %>% arrange(group, particle, segmenttronos)
    durations <- df_cluster$frame.e - df_cluster$frame.b + 1
    ret <- rep(df_cluster$cluster, durations)
    return(data.frame(df_segment, cluster = ret))
}

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

calculate.spectrum <- function(y, time = NULL) {
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

calculate.freqs <- function(y, time = NULL, sampling.time = 1) {
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

calculate.freqs.max <- function(spx, spy) {
    return(spx[spy == max(spy, na.rm = TRUE)][1])
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

spectrum.slope <- function(y) {
    wave.out <- morlet(y, 1:length(y), p2 = 8, dj = 0.1, siglvl = 0.999)
    # TODO: change the frequency scanning!
    prd <- wave.out$period[apply(Re(wave.out$wave)[,1:50], 1, which.max)]
    mod <- lm(prd ~ c(1:length(y)))
    mod.slope <- summary(mod)[[4]][2]
    return(c(mod.slope, mod.slope / abs(mod.slope)))
}

segment.timeseries <-
    function(df,
             coords = c("X", "Y"),
             sampling.time,
             updateProgress = NULL) {

        df.proc <- df %>%
            group_by(group, particle) %>%
            mutate(group = factor(group), particle = factor(particle))

        # TODO: segment only one variable
        df.segm <- df.proc %>%
            mutate(segmenttronos = segment.lavielle(!!sym(coords[1]), !!sym(coords[2])),
                frame.init = first(frame),
                frame.seg = c(1:length(!!sym(coords[1])))
            )

        if (is.function(updateProgress)) {
            updateProgress("Series segmented")
        }

        df.spec <- df.segm %>%
            group_by(group, particle, segmenttronos) %>%
            mutate(spec.s = calculate.spectrum(y),
                   spec.f = calculate.freqs(y, sampling.time = sampling.time)
            ) %>% top_n(1, spec.s) %>% ungroup() %>% arrange(group, particle, segmenttronos) %>%
             select(spec.s, spec.f)

        if (is.function(updateProgress)) {
            updateProgress("Spectrums calculated")
        }

        if (length(coords) == 1) {
            df.summ <- df.segm %>%
            group_by(group, particle, segmenttronos) %>%
            arrange(group, particle, segmenttronos) %>%
            summarise(
                mean.c1 = mean(!!sym(coords[1]), na.rm = TRUE),
                sd.c1 = sd(!!sym(coords[1]), na.rm = TRUE),
                frame.b = first(frame.seg),
                frame.e = last(frame.seg),
                begin = first(frame.seg) + first(frame.init) - 1,
                end = last(frame.seg) + first(frame.init) - 1
            ) %>% arrange(group, particle, segmenttronos)
        } else {
            df.spec.ms.m <- df.segm %>%
                group_by(group, particle, segmenttronos) %>%
                bow(tie(spec.m, spec.ms) := spectrum.slope(y)) %>%
                ungroup() %>% arrange(group, particle, segmenttronos) %>% select(spec.m, spec.ms)
            df.summ <- df.segm %>%
                group_by(group, particle, segmenttronos) %>%
                summarise(
                    mean.c1 = mean(!!sym(coords[1]), na.rm = TRUE),
                    mean.c2 = mean(!!sym(coords[2]), na.rm = TRUE),
                    sd.c1 = sd(!!sym(coords[1]), na.rm = TRUE),
                    sd.c2 = sd(!!sym(coords[2]), na.rm = TRUE),
                    frame.b = first(frame.seg),
                    frame.e = last(frame.seg),
                    begin = first(frame.seg) + first(frame.init) - 1,
                    end = last(frame.seg) + first(frame.init) - 1
                ) %>% arrange(group, particle, segmenttronos)
            df.summ <- cbind(df.summ, df.spec, df.spec.ms.m)
        }

        if (is.function(updateProgress)) {
            updateProgress("Summary computed")
        }

        df.ret <- df.summ %>% ungroup()
        return(list(as.data.frame(df.ret), as.data.frame(df.segm)))
    }

segment.timeseries.multiv <-
    function(df,
             coords,
             sampling.time,
             updateProgress = NULL) {
        df.ret <- data.frame()
        df.segm <- data.frame()
        segments <- c()
        groups <- df %>% select(group) %>% distinct()
        for (g in groups[, 1]) {
            particles <-
                df %>% filter(group == g) %>% select(particle) %>% distinct()
            for (p in particles[, 1]) {
                # Parsing time series
                seq.temp <-
                    df %>% filter(group == g) %>% filter(particle == p)

                time.zero <- seq.temp[["frame"]][1]
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
                curr.seg <- data.frame(segments.calc$segments)
                nsegments <- nrow(segments.calc$segments)
                names(curr.seg)[2] <- "begin"
                time.vec <- as.vector(t(curr.seg[,c(2,3)]))
                time.vec[seq(length(time.vec)) %% 2 == 1] <- time.vec[seq(length(time.vec)) %% 2 == 1] - 1
                time.vec <- c(time.vec, nrow(seq.temp.filtered))
                segments <- c(segments, rep(c(rbind(1:nsegments, NA)), c(diff(time.vec))))
                curr.seg <- curr.seg %>% select(begin, end)
                curr.seg$group <- factor(g)
                curr.seg$state <- factor(c(1:nsegments))
                curr.seg$segmenttronos <- factor(c(1:nsegments))
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
                    
                    spec.s = calculate.spectrum(seg.temp[[coords[1]]], seg.temp[["frame"]])
                    spec.f = calculate.freqs(seg.temp[[coords[1]]], seg.temp[["frame"]], sampling.time)
                    
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
                df.ret <- rbind(df.ret, curr.seg)
                if (is.function(updateProgress)) {
                    text <- paste0("Particle: ", p, ", Group: ", g)
                    updateProgress(text)
                }
            }
        }
        df.segm <- df %>% mutate(segmenttronos = segments) %>% na.omit()
        return(list(as.data.frame(df.ret), as.data.frame(df.segm)))
    }

clusterize.segments <-
    function(class.diff, grouping, particle, segments, updateProgress = NULL) {
        if (grouping != "none" && particle != "none") {
            drop.cols <-
                c("segmenttronos",
                  "begin",
                  "end",
                  "group",
                  "particle",
                  "state")
            
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

null.if.empty <- function(x) {
    if(length(x) == 0) {
        return(NULL)
    }

    return(x)
}

motifs.per.segment <- function(df.segm, var, nmain, nquery, window, pct, thr) {
    segments.motifs.list <- df.segm %>% filter(group == nmain) %>%
            group_by(cluster) %>% dplyr::filter(dplyr::n() * 0.5 >= window) %>%
            group_map( ~ motif.discovery(
                .[[var]],
                win = window,
                null.if.empty(df.segm[df.segm$group == nquery && df.segm$cluster == .$cluster, var]),
                pct = pct,
                thr = thr
            ))

    return(segments.motifs.list)
}

find.pc.causality <- function(X, pval, updateProgress = NULL) {
  if (is.function(updateProgress)) {
    updateProgress("Step 1/3: Computing PC algorithm")
  }
  stuffStat <- list(C = cor(as.matrix(X)), n = nrow(as.matrix(X)))
  pc.fit <- pc(stuffStat, indepTest = gaussCItest, p = ncol(as.matrix(X)), alpha = pval)
  as(pc.fit, "matrix")
}

find.vlgc.causality <- function(X, nlags, pval, updateProgress = NULL) {
  if (is.function(updateProgress)) {
    updateProgress("Step 3/3: Computing VL-GC algorithm")
  }
  out <- multipleVLGrangerFunc(as.matrix(X), maxLag = as.numeric(nlags), alpha = pval)
  out$adjMat
}

find.vlte.causality <- function(X, nlags, pval, updateProgress = NULL) {
  if (is.function(updateProgress)) {
    updateProgress("Step 2/3: Computing VL-TE algorithm")
  }
  out <- multipleVLTransferEntropy(as.matrix(X), maxLag = as.numeric(nlags), autoLagflag = TRUE, alpha = pval)
  out$adjMat
}

causality.global <- function(df, vars, pval, nlags, group, progress = NULL){
    # Mutate variables of signal and morphology

    df_filter <- df[,c("particle", "group", vars)] %>% filter(group == group) %>% select(-group)
    nparts <- nrow(df_filter %>% dplyr::select(particle) %>% distinct())
    updateProgress <- function(detail = NULL) {
      progress$inc(amount = 1/(nparts*3), detail = detail)
    }

    # Case 1: PC causality
    list.causality.pc <- df_filter %>% dplyr::group_by(particle) %>%
        group_map( ~ find.pc.causality(., pval, updateProgress))

    adj.mat.pc <- round(Reduce('+', list.causality.pc) / nparts, 2)

    # Case 2: VLTE causality
    list.causality.vlte <- df_filter %>% dplyr::group_by(particle) %>%
        group_map( ~ find.vlte.causality(., nlags, pval, updateProgress))

    adj.mat.vlte <- round(Reduce('+', list.causality.vlte) / nparts, 2)

    # Case 3: VLGC causality
    list.causality.vlgc <- df_filter %>% dplyr::group_by(particle) %>%
        group_map( ~ find.vlgc.causality(., nlags, pval, updateProgress))

    adj.mat.vlgc <- round(Reduce('+', list.causality.vlgc) / nparts, 2)

    return(list(adj.mat.pc, adj.mat.vlte, adj.mat.vlgc))
}

causality.plot <- function(adj.mat, threshold, size, names){
    adj.mat <- apply(adj.mat, 1, function(x) ifelse(x >= threshold, x, 0))
    out.network <- as.network(adj.mat,mode="directed", weighted = T, ignore.eval = FALSE, names.eval = "weights", directed = T)
    network.vertex.names(out.network) <- names
    tryCatch({
        ggnet2(out.network, size = size, label = TRUE, label.size = 5, arrow.size = 12, arrow.gap = 0.05, edge.size = "weights", edge.label = "weights")
    },error=function(cond){
        ggnet2(out.network, size = size, label = TRUE, label.size = 5, arrow.size = 12, arrow.gap = 0.05)
    })
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
    drop.cols <- c("subsample_ind", "subsample_ind.y", "state", "group", "particle", "segmenttronos")
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

    partial_bundle(toWebGL(subplot(ggplotly(p1, height = 500), ggplotly(p2, height = 500), widths = c(0.7, 0.3), titleX=T)))
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
             coords,
             sampling.time,
             updateProgress = NULL) {

        df_freqs <- df %>% group_by(particle) %>%
            mutate(spec.s = calculate.spectrum(!!sym(coords[1])),
                   spec.f = calculate.freqs(!!sym(coords[1]),sampling.time = sampling.time)) %>%
            na.omit() %>% mutate(spec.s = normalize.z(spec.s))

        return(as.data.frame(df_freqs))
    }

spectrogram.calculate <- function(y, frame) {
    # TODO: speedup the function (not really that much is possible)
    wave.out <- morlet(y, 1:length(y), p2 = 8, dj = 0.1)
    Power <- wave.out$Power
    # TODO: Select the interval for the morlet spectrogram
    Signif <- as.data.frame(Power[,1:50])
    names(Signif) <- wave.out$period[1:50]

    Signif <- data.frame(t(apply(Signif, 1, function(x) head(names(Signif)[order(-x)],10))))
    Signif <- cbind(Signif, data.frame(frame = frame))
    Signif <- Signif %>% pivot_longer(!frame, names_to = "Y", values_to = "spec") %>%
                mutate_each_(funs(as.numeric), "spec")
    return(as.data.frame(Signif))
}

spectrogram.global <-
    function(df.segm,
             coords,
             updateProgress = NULL) {

        df.spec <- df.segm %>%
            group_by(group, particle) %>% mutate(SVar = !!sym(coords[1])) %>%
            do(spectrogram.calculate(.$SVar, .$frame)) %>%
            na.omit() %>% mutate(label = group)

        return(df.spec)
    }

spectrogram.per.segment <-
    function(df.segm,
             coords,
             updateProgress = NULL) {
        df.spec <- df.segm %>%
            group_by(group, particle, cluster) %>% mutate(SVar = !!sym(coords[1])) %>%
            do(possibly(spectrogram.calculate, otherwise = data.frame(frame = NA, Y = NA, spec = NA))(.$SVar, .$frame)) %>%
            na.omit() %>% mutate(label = group)
        return(df.spec)
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
             sampling.time,
             freq = FALSE,
             individual = FALSE) {

        if (freq) {
            p1 <-
                ggplot(df_freqs, aes(spec.f, spec.s), colour = group) +
                xlab("Frequency")
        } else {
            p1 <-
                ggplot(df_freqs, aes((1 / spec.f) / 60, spec.s), colour = group) +
                xlab("Period (minutes)")
        }
        p1 <- p1 + stat_smooth(
                method = "loess",
                span = 0.1,
                se = TRUE,
                aes(
                    color = group,
                    fill = group
                ),
                alpha = 0.3
            ) +
            theme_bw() +
            ylab("Normalized spectral density") +
            scale_x_continuous(trans = "log10")

        plt <- ggplotly(p1, height = 500, width = 600)

        if (individual) {
            npar <- df_freqs %>% ungroup() %>% select(particle) %>% distinct()
            p1 <- p1 + facet_wrap(particle ~ .) +
            theme(panel.spacing.x = unit(1, "lines"), panel.spacing.y = unit(2, "lines"))
            plt <- ggplotly(p1, height = 75 * nrow(npar), width = 600)
        }

        return(partial_bundle(toWebGL(plt)))
    }

spectral.plot.segments <-
    function(df_seg_freqs,
             sampling.time,
             freq = FALSE) {

        nclust <- df_seg_freqs %>% ungroup() %>% select(cluster) %>% distinct()

        if(is.null(df_seg_freqs)) {
            return(NULL)
        }

        if (freq) {
            p1 <-
                ggplot(df_seg_freqs, aes(spec.f, spec.s), colour = group) +
                xlab("Frequency")
        } else {
            p1 <-
                ggplot(df_seg_freqs, aes((1 / spec.f) / sampling.time, spec.s), colour = group) +
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
            theme(panel.spacing.x = unit(1, "lines"), panel.spacing.y = unit(2, "lines"))

        plt <- ggplotly(p1, height = 75 * nrow(nclust), width = 600)

        return(partial_bundle(toWebGL(plt)))
    }

heatmap.plot <-
    function(df_specgram,
             freq = FALSE) {

        # x <- list(
        #     title = "Time"
        # )
        # y <- list(
        #     title = "Period (min)"
        # )
        # df <- data.frame(df_specgram)
        # print(df)

        # plt <- df %>%
        #         split(df$group) %>%
        #         purrr::map(~{
        #             plot_ly(data = .x, x = .x$frame, y = .x$spec, name = .x$group, colors = colorRamp(specCols)) %>%
        #             add_histogram2dcontour(ncontours=40, contours = list(coloring='heatmap'), line = list(width = 0)) %>%
        #             layout(yaxis = y, xaxis = x, annotations = list(
        #                                                         text = .x$group[1],
        #                                                         xref = "paper",
        #                                                         yref = "paper",
        #                                                         yanchor = "bottom",
        #                                                         xanchor = "center",
        #                                                         align = "center",
        #                                                         x = 0.5,
        #                                                         y = 1,
        #                                                         showarrow = FALSE
        #                                                     ),
        #                     plot_bgcolor='#5E4FA2')
        #         }) %>% subplot(nrows = length(.), margin = 0.05, shareX = TRUE)

        if (freq) {
            p1 <-
                df_specgram %>% ggplot(aes(x=frame, y=(1/spec)/60) ) +
                xlab("Frequency")
        } else {
            p1 <-
                df_specgram %>% ggplot(aes(x=frame, y=spec)) +
                xlab("Period")
        }
            p1 <- p1 + stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
                scale_fill_gradientn(colors = specCols) +
                facet_grid(label ~ .) +
                ylab("Period") +
                theme_bw()

        return(p1)
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
            p1 <- p1 + stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
                facet_grid(cluster ~ label * .) +
                scale_fill_gradientn(colors = specCols) +
                xlab("Time") +
                theme_bw()

        plt <- p1

        return(plt)
    }

heatmap.veloc.plot <-
    function(df_vel) {
        # x <- list(
        #     title = "Time"
        # )
        # y <- list(
        #     title = "Velocity"
        # )

        # # steps_x <- list()
        # # for (i in 1:50) {
        # #     step <- list(args = list("xbins.size", i/10),
        # #                     label = i/10,
        # #                     method = "restyle",
        # #                     value = i/10
        # #                     )

        # #     steps_x[[i]] = step
        # # }

        # # steps_y <- list()
        # # for (i in 1:50) {
        # #     step <- list(args = list("ybins.size", i/10),
        # #                     label = i/10,
        # #                     method = "restyle",
        # #                     value = i/10
        # #                     )

        # #     steps_y[[i]] = step
        # # }

        # df <- data.frame(df_vel)

        # plt <- df %>%
        #         split(df$group) %>%
        #         purrr::map(~{
        #             plot_ly(data = .x, x = .x$frame, y = .x$vel.ma, name = .x$group, colors = colorRamp(specCols)) %>%
        #             add_histogram2dcontour(contours = list(coloring='heatmap'), line = list(width = 0)) %>%
        #             layout(yaxis = y, xaxis = x, annotations = list(
        #                                                         text = .x$group[1],
        #                                                         xref = "paper",
        #                                                         yref = "paper",
        #                                                         yanchor = "bottom",
        #                                                         xanchor = "center",
        #                                                         align = "center",
        #                                                         x = 0.5,
        #                                                         y = 1,
        #                                                         showarrow = FALSE
        #                                                     ),
        #                     plot_bgcolor='#5E4FA2'
        #                     # sliders = list(
        #                     #             list(
        #                     #                 active = 2,
        #                     #                 currentvalue = list(prefix = "X Bin size: "),
        #                     #                 pad = list(t = 20),
        #                     #                 steps = steps_x),
        #                     #             list(
        #                     #                 active = 2,
        #                     #                 currentvalue = list(prefix = "Y Bin size: "),
        #                     #                 pad = list(t = 100),
        #                     #                 steps = steps_y)
        #                     #             )
        #                     # )
        #             )
        #         }) %>% subplot(nrows = length(.), margin = 0.05, shareX = TRUE)

        names(df_vel)[names(df_vel) == "group"] <- 'label.tronos'

        p1 <- df_vel %>% ggplot(aes(x= frame, y= vel.ma)) +
                xlab("Time") +
                stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
                scale_fill_gradientn(colors = specCols) +
                facet_grid(label.tronos ~ .) +
                ylab("Velocity") +
                theme_bw()

        return(p1)
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
            stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
            facet_grid(cluster ~ group * .) +
            scale_fill_gradientn(colors = specCols) +
            xlab("Time") +
            theme_bw()

        plt <- p1

        return(plt)
    }

spectrum.significance <-
    function(df_freqs,
             range) {
        df_freqs <- df_freqs %>% mutate(spec.f = spec.f / max(spec.f)) %>%
            filter(spec.f < range[2], spec.f > range[1])
        
        m1.r <-
            glm(factor(group) ~ spec.f * spec.s ,
                data = df_freqs,
                family = "binomial")
        m2.r <-
            glm(factor(group) ~ spec.f + spec.s ,
                data = df_freqs,
                family = "binomial")
        m3.r <-
            glm(factor(group) ~ 1 ,
                data = df_freqs,
                family = "binomial")
        anova(m1.r, m2.r, m3.r, test = "Chisq")
    }

spectrum.per.segment <-
    function(df.segm,
             coord,
             sampling.time,
             updateProgress = NULL) {

        if (is.function(updateProgress)) {
            updateProgress("Preparing data")
        }

        df.spec <- data.frame(df.segm) %>% group_by(group, particle) %>% dplyr::filter(dplyr::n() >= 10) %>%
            group_by(group, particle, cluster) %>% dplyr::filter(dplyr::n() >= 5) %>%
            mutate(spec.s = calculate.spectrum(y),
                   spec.f = calculate.freqs(y, sampling.time = sampling.time)
            ) %>% na.omit()

        if (is.function(updateProgress)) {
            updateProgress("Spectrum computed")
        }

        return(df.spec)
    }

velocity.per.segment <-
    function(df.segm,
             coords,
             vel_ma,
             updateProgress = NULL) {

        if (is.function(updateProgress)) {
            updateProgress("Preparing data")
        }


        df.vels <- df.segm %>% group_by(group, particle) %>% dplyr::filter(dplyr::n() >= 10) %>%
            group_by(group, particle, cluster) %>% dplyr::filter(dplyr::n() >= vel_ma * 2) %>%
            mutate(vel.ma = calculate.velocity.ma(!!sym(coords[1]), !!sym(coords[2]), vel_ma)) %>%
            na.omit() %>% mutate(label = group)

        if (is.function(updateProgress)) {
            updateProgress("Velocities computed")
        }

        vels.segs.grouped <- df.vels %>% ungroup()
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

calc.velocities <- function(df, coords, vel_ma) {
    df_vt.smooth <- df %>% group_by(group, particle) %>%
        mutate(vel.ma = calculate.velocity.ma(!!sym(coords[1]),!!sym(coords[2]), vel_ma)) %>%
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

significant.velocities <- function(df_vt.smooth) {
    wilcox.test(df_vt.smooth$vel.ma ~ df_vt.smooth$group, conf.int = TRUE)
}

plot.velocities <-
    function(df_vt.smooth, individual = FALSE) {
        p1 <-
            ggplot(df_vt.smooth, aes(
                y = vel.ma,
                fill = group,
                x = group
            )) +
            geom_violin(trim = TRUE, alpha = 0.5) +
            geom_boxplot(width = 0.1,
                         trim = TRUE,
                         outlier.shape = NA) +
            theme_bw() +
            ylab("Velocity")

        if (individual) {
            p1 <- p1 + facet_wrap(particle ~ .)
        }

        plt <- ggplotly(p1, autosize = TRUE, height = 500, width = 600)
        # plt <- ggarrange(p1, p2, nrow = 2)
        return(plt)
    }

plot.velocities.segment <-
    function(df_seg_vels) {

        nclust <- df_seg_vels %>% ungroup() %>% select(cluster) %>% distinct()

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

        plt <- ggplotly(p2, height = 75 * nrow(nclust), width = 600)
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

        plt <- ggplotly(p1, height = 600, width = 600)
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
        return(partial_bundle(toWebGL(plt)))
    }

clustering.variables <- function(df, variables) {
    # Clustering to retrieve morphology classes
    classes <- Mclust(df %>% select(variables))
    return(classes)
}

clustering.plot <- function(classes) {
    return(plot.Mclust(classes))
}