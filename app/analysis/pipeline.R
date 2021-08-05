# ChroMo
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

create.df.clusters <- function(df_cluster, df_segment) {
    df_cluster <- df_cluster %>% arrange(group, particle, segmentchromo)
    df_segment <- df_segment %>% arrange(group, particle, segmentchromo)
    durations <- df_cluster$frame.e - df_cluster$frame.b + 1
    ret <- rep(df_cluster$cluster, durations)
    return(data.frame(df_segment, cluster = ret))
}

segment.timeseries <-
    function(df,
             coords = c("X", "Y"),
             sampling.time,
             updateProgress = NULL) {

        if (length(coords) != 2) {
            coords <- c(coords[1], coords[1])
        }
        df.proc <- df %>%
            group_by(group, particle) %>%
            mutate(group = factor(group), particle = factor(particle))

        # TODO: segment only one variable
        df.segm <- df.proc %>%
            mutate(segmentchromo = segment.lavielle(!!sym(coords[1]), !!sym(coords[2])),
                frame.init = first(frame),
                frame.seg = c(1:length(!!sym(coords[1])))
            )

        if (is.function(updateProgress)) {
            updateProgress("Series segmented")
        }

        df.spec <- df.segm %>%
            group_by(group, particle, segmentchromo) %>%
            mutate(spec.s = calculate.spectrum(y),
                   spec.f = calculate.freqs(y, sampling.time = sampling.time)
            ) %>% top_n(1, spec.s) %>% ungroup() %>% arrange(group, particle, segmentchromo) %>%
             select(spec.s, spec.f)

        if (is.function(updateProgress)) {
            updateProgress("Spectrums calculated")
        }

        df.spec.ms.m <- df.segm %>%
            group_by(group, particle, segmentchromo) %>%
            bow(tie(spec.m, spec.ms) := spectrum.slope(y)) %>%
            ungroup() %>% arrange(group, particle, segmentchromo) %>% select(spec.m, spec.ms)
        df.summ <- df.segm %>%
            group_by(group, particle, segmentchromo) %>%
            summarise(across(coords,
                    list(mean=~mean(.x, na.rm=TRUE),
                         sd=~sd(.x, na.rm=TRUE),
                         median=~median(.x, na.rm=TRUE))
                ),
                frame.b = first(frame.seg),
                frame.e = last(frame.seg),
                begin = first(frame.seg) + first(frame.init) - 1,
                end = last(frame.seg) + first(frame.init) - 1
            ) %>% arrange(group, particle, segmentchromo)
        df.summ <- cbind(df.summ, df.spec, df.spec.ms.m)

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
                curr.seg$segmentchromo <- factor(c(1:nsegments))
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
                        seq.temp.filtered[beg:end, ] %>% summarise_all(list(
                            mean=~mean(.x, na.rm = TRUE),
                            sd=~sd(.x, na.rm = TRUE),
                            median=~median(.x, na.rm = TRUE)
                        ))

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
        df.segm <- df %>% mutate(segmentchromo = segments) %>% na.omit()
        return(list(as.data.frame(df.ret), as.data.frame(df.segm)))
    }

clusterize.segments <-
    function(class.diff, grouping, particle, segments, updateProgress = NULL) {
        if (grouping != "none" && particle != "none") {
            drop.cols <-
                c("segmentchromo",
                  "begin",
                  "end",
                  "frame.e",
                  "frame.b",
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


motifs.per.segment <- function(df.segm, var, nmain, nquery, window, pct, thr, motif_amount = 3L) {
    segments.motifs.list <- df.segm %>% filter(group %in% nmain) %>%
            group_by(cluster) %>% dplyr::filter(dplyr::n() * 0.5 >= window) %>%
            group_map( ~ motif.discovery(
                .[[var]],
                win = window,
                NULL,
                pct = pct,
                thr = thr,
                time = .[["frame"]],
                group = .[["group"]],
                nmotifs = motif_amount
            ))

    return(segments.motifs.list)
}

segmentation.summary <- function(segments) {
    drop.cols <- c("subsample_ind", "subsample_ind.y", "state", "group", "particle", "segmentchromo")
    segments <- segments %>% select(-one_of(drop.cols))
    return(segments %>% group_by(cluster) %>% summarise_all(mean))
}

segmentation.grouped.summary <- function(segments) {
    drop.cols <- c("subsample_ind", "subsample_ind.y", "state", "particle", "segmentchromo")
    segments <- segments %>% select(-one_of(drop.cols))
    return(segments %>% group_by(group, cluster) %>% summarise_all(mean))
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
            mutate(vel.ma = calculate.velocity.ma(., coords, vel_ma)) %>%
            na.omit() %>% mutate(label = group)

        if (is.function(updateProgress)) {
            updateProgress("Velocities computed")
        }

        vels.segs.grouped <- df.vels %>% ungroup()
        return(vels.segs.grouped)
    }

density.per.segment <-
    function(df.segm,
             coords,
             vel_ma,
             updateProgress = NULL) {

        if (is.function(updateProgress)) {
            updateProgress("Preparing data")
        }


        df.vels <- df.segm %>% group_by(group, particle) %>% dplyr::filter(dplyr::n() >= 10) %>%
            group_by(group, particle, cluster) %>% dplyr::filter(dplyr::n() >= vel_ma * 2) %>%
            mutate(vel.ma = calculate.ma(!!sym(coords[1]), vel_ma)) %>%
            na.omit() %>% mutate(label = group)

        if (is.function(updateProgress)) {
            updateProgress("Densities computed")
        }

        vels.segs.grouped <- df.vels %>% ungroup()
        return(vels.segs.grouped)
    }


motif.discovery <-
    function(df,
             win = 12,
             q = NULL,
             pct = 1,
             thr = 0.98,
             time = NULL,
             group = NULL,
             nmotifs = 3L) {

        mp <- tsmp::compute(df,
            windows = win,
            query = NULL,
            sample_pct = pct,
            threshold = thr,
            n_jobs = 1L)
        motifs <- tsmp::motifs(mp, k = nmotifs, neighbor_count = 100L)
        discords <- tsmp::discords(mp, k = nmotifs, neighbor_count = 100L)

        motifs.idx <- motifs[["motif"]][["motif_idx"]]
        discord.idx <- discords[["discord"]][["discord_idx"]]
        motifs.neigh <- motifs[["motif"]][["motif_neighbor"]]
        discord.neigh <- discords[["discord"]][["discord_neighbor"]]
        motifs.loc <- lapply(seq_along(motifs.idx), function(i) c(motifs.idx[[i]], motifs.neigh[[i]]))
        discord.loc <- lapply(seq_along(discord.idx), function(i) c(discord.idx[[i]], discord.neigh[[i]]))
        motifs.group <- lapply(motifs.loc, function(i) group[i])
        discord.group <- lapply(discord.loc, function(i) group[i])
        motifs.loc <- lapply(motifs.loc, function(i) time[i])
        discord.loc <- lapply(discord.loc, function(i) time[i])

        motifs.locations <- lapply(
            seq_along(motifs.loc),
            function(i) rbind(
                data.frame(motif = i, location = motifs.loc[i], group = motifs.group[i])) %>%
                    dplyr::rename(location = 2, group = 3)
            ) %>% bind_rows()


        discord.locations <- lapply(
            seq_along(discord.loc),
            function(i) rbind(
                data.frame(motif = i, location = discord.loc[i], group = discord.group[i])) %>%
                    dplyr::rename(location = 2, group = 3)
            ) %>% bind_rows()

        return(list(motifs, motifs.locations, discord.locations))
    }

clustering.variables <- function(df, variables) {
    # Clustering to retrieve morphology classes
    classes <- Mclust(df %>% select(variables))
    return(classes)
}

reapply.clustering <- function(df, p.corr = NULL) {
    print(df)
    df.output <- df
    drop.cols <-
                c("segmentchromo",
                  "frame.init",
                  "frame.seg",
                  "begin",
                  "end",
                  "frame.e",
                  "frame.b",
                  "particle")

    class.diff.cluster <-
                df %>% select(-one_of(drop.cols))

    cluster.no <- df %>% select(cluster) %>% distinct()
    group.no <- df %>% select(group) %>% distinct()
    cluster.changes <- df %>% select(group, cluster) %>% distinct()
    cluster.changes <- cbind(cluster.changes, data.frame(change = cluster.changes$cluster))
    for (i in group.no[, 1]) for (j in group.no[, 1]) {
        if (i == j) next
        for (c in cluster.no) {
            current.cluster <- class.diff.cluster %>% filter(cluster == c)
            current.groups <- current.cluster %>% filter(group %in% c(i, j))
            print(current.groups)
            mod.1 <-
                glm(group ~ . - cluster,
                    data = current.groups,
                    family = "binomial")
            mod.2 <-
                glm(group ~ 1,
                    data = current.groups,
                    family = "binomial")

            test <- anova(mod.1, mod.2)
            print(test)

            # p.val <- summary(test)[[1]][["Pr(>F)"]]

            # if (p.val < 0.05) {
            #     print("Significant segment differences found!")
            #     max.clust <- max(cluster.changes$change)
            #     df.output <- df.output %>% mutate(cluster = ifelse(group == j, cluster + max.clust, cluster))
            #     cluster.changes <- rbind(cluster.changes, data.frame(group = j, cluster = c, change = c + n.clusters))
            # } else {
            #     print("No significant differences found!")
            #     go.back <- cluster.changes %>% filter(group == i && cluster == c) %>% select(change)
            #     df.output <- df.output %>% mutate(cluster = ifelse(group == j, go.back, cluster))
            # }
        }
    }

    return(df.output)
}