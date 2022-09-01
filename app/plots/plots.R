segmentation.plot <- function(segments, grouping = "none") {
    segno <- segments %>% group_by(group) %>% count(cluster) %>% mutate(freq = n / sum(n))

    p1 <-
        ggplot(data = segments, aes(x = begin, y = cluster, color = group)) +
        facet_grid(group ~ .) +
        geom_segment(aes(x = begin, xend = end, yend = cluster), alpha = 0.2) +
        geom_point(aes(x = begin), size = 3, alpha = 0.2) +
        geom_point(aes(x = end), size = 3, alpha = 0.2) +
        theme_bw() +
        xlab("Time") +
        ylab("Cluster") +
        theme(legend.position = "none")

    p2 <- ggplot(data = segno, aes(x = cluster, y = freq, fill = group)) +
        facet_grid(group ~ .) +
        geom_col(alpha = 0.5) +
        coord_flip() +
        theme_bw() +
        xlab("Cluster") +
        ylab("Count") +
        theme(legend.position = "none")

    ggarrange(p1, p2, widths = c(0.7, 0.3))
}

segmentation.plotly <- function(segments, grouping = "none") {
    segno <- segments %>% group_by(group) %>% count(cluster) %>% mutate(freq = n / sum(n))

    p1 <-
        ggplot(data = segments, aes(x = begin, y = cluster, color = group)) +
        facet_grid(group ~ .) +
        geom_segment(aes(x = begin, xend = end, yend = cluster), alpha = 0.2) +
        geom_point(aes(x = begin), size = 3, alpha = 0.2) +
        geom_point(aes(x = end), size = 3, alpha = 0.2) +
        theme_bw() +
        xlab("Time") +
        ylab("Cluster") +
        theme(legend.position = "none")

    p2 <- ggplot(data = segno, aes(x = cluster, y = freq, fill = group)) +
        facet_grid(group ~ .) +
        geom_col(alpha = 0.5) +
        coord_flip() +
        theme_bw() +
        xlab("Cluster") +
        ylab("Count") +
        theme(legend.position = "none")

    partial_bundle(toWebGL(subplot(ggplotly(p1, height = 500), ggplotly(p2, height = 500), widths = c(0.7, 0.3), titleX=T)))
}

morlet.plot <- function(df, plot.var, sig = 0.95, ylab = "Y coordinate") {
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
        useRaster = TRUE,
        crn.lab = ylab
    )
}

spectral.plot <-
    function(df_freqs,
             sampling.time,
             freq = FALSE,
             individual = FALSE, 
             x_ticks = 10) {

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
            ylab("Normalized spectral density")

        x_lim <- layer_scales(p1)$x$get_limits()
        p1 <- p1 + scale_x_continuous(trans = "log10", 
                                      breaks = 10^seq(log10(x_lim[1]), 
                                                            log10(x_lim[2]), 
                                                            length.out = x_ticks))

        plt <- ggplotly(p1, height = 500, width = 600)

        if (individual) {
            npar <- df_freqs %>% ungroup() %>% select(particle) %>% distinct()
            p1 <- p1 + facet_wrap(particle ~ .) +
            theme(panel.spacing.x = unit(1, "lines"), panel.spacing.y = unit(2, "lines"))
            plt <- ggplotly(p1, height = 75 * nrow(npar), width = 600)
        }

        return(p1)
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

        return(p1)
    }

heatmap.plot <-
    function(df_specgram,
             freq = FALSE) {

        if (freq) {
            p1 <-
                df_specgram %>% ggplot(aes(x=frame, y=(1/spec)/60) ) +
                ylab("Frequency")
        } else {
            p1 <-
                df_specgram %>% ggplot(aes(x=frame, y=spec)) +
                ylab("Period")
        }
            p1 <- p1 + stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
                scale_fill_gradientn(colors = specCols) +
                facet_grid(label ~ .) +
                xlab("Time") +
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
    function(df_vel, colname="Velocity") {
        names(df_vel)[names(df_vel) == "group"] <- 'label.chromo'

        p1 <- df_vel %>% ggplot(aes(x= frame, y= vel.ma)) +
                xlab("Time") +
                stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
                scale_fill_gradientn(colors = specCols) +
                facet_grid(label.chromo ~ .) +
                ylab(colname) +
                theme_bw()

        return(p1)
    }

heatmap.veloc.segmented.plot <-
    function(df_vel_segmented,
             freq = FALSE,
             colname="Velocity") {

        if(is.null(df_vel_segmented)) {
            return(NULL)
        }

        p1 <-
            df_vel_segmented %>% ggplot(aes(x=frame, y=vel.ma )) +
            ylab(colname) +
            stat_density_2d(aes(fill = ..ndensity..), geom = "raster", contour = FALSE) +
            facet_grid(cluster ~ group * .) +
            scale_fill_gradientn(colors = specCols) +
            xlab("Time") +
            theme_bw()

        plt <- p1

        return(plt)
    }

plot.velocities <-
    function(df_vt.smooth, individual = FALSE, colname="Velocity") {
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
            ylab(colname)

        if (individual) {
            p1 <- p1 + facet_wrap(particle ~ .)
        }

        plt <- ggplotly(p1, autosize = TRUE, height = 500, width = 600)
        return(p1)
    }

plot.velocities.segment <-
    function(df_seg_vels, colname="Velocity") {

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
            ylab(colname) +
            theme(panel.spacing.x = unit(1, "lines"), panel.spacing.y = unit(1, "lines"))

        plt <- ggplotly(p2, height = 100 * nrow(nclust), width = 600)
        return(p2)
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
        return(p1)
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
            facet_wrap(cluster ~ ., scale = "free") +
            theme_bw() +
            ylab("Mean Squared Displacement") +
            xlab("Lag")

        plt <- ggplotly(p1, height = 500, width = 600)
        return(p1)
    }

locations.motifs.plot <- function(locations) {

    nlocations <- nrow(locations %>% select(motif) %>% distinct())

    p1 <- locations %>% ggplot(aes(x = group, y = location, color = group, label = particle)) +
            geom_violin() +
            geom_jitter() +
            facet_wrap(motif ~ .) +
            xlab("") +
            ylab("Time location") +
            theme_bw() +
            theme(legend.position = "none")

    plt <- ggplotly(p1, height = 200, width = 200 * nlocations)

    return(p1)
}

causality.plot <- function(adj.mat, threshold, names){
    adj.mat <- t(apply(adj.mat, 1, function(x) ifelse(x >= threshold, x, 0)))
    net_build <- igraph_to_networkD3(graph_from_adjacency_matrix(adj.mat, mode = "directed", weighted = TRUE))
    net_build$nodes$name <- names
    net_build$nodes$group <- 1
    net_build$links$value <- net_build$links$value * 2
    tryCatch({
        forceNetwork(Links = net_build$links, Nodes = net_build$nodes,
            Source = "source",Target = "target", Value = "value", NodeID = "name",
            Group = "group", opacity = 0.8, arrows = TRUE, fontSize = 12, fontFamily = "Arial",
            opacityNoHover = 1, zoom = TRUE)
    },error=function(cond){
        forceNetwork(Links = data.frame(source = 0, target = 1, value = 0),
            Nodes = net_build$nodes,
            Source = "source",Target = "target", Value = "value", NodeID = "name",
            Group = "group", opacity = 0.8, arrows = TRUE, fontSize = 12, fontFamily = "Arial",
            opacityNoHover = 1, zoom = TRUE)
    })
}

clustering.plot <- function(classes) {
    return(plot.Mclust(classes))
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