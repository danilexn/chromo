example.data <- read.csv("../examples/example.csv")
example.data.huge <- read.csv("../examples/example_huge.csv")
vals <- reactiveValues(count = 0)

server <- function(input, output, session) {
  isolate(vals$count <- vals$count + 1)
  fileUploaded <- FALSE
  setBookmarkExclude(c("run_segmentation", "run_msd", "run_motifs", "run_causal", "run_spectrum", "run_velocities"))
  restore_segment_variables <- reactiveVal()
  restore_df_vars_motifs <- reactiveVal()
  restore_df_query_causality <- reactiveVal()
  restore_df_sequence_motifs <- reactiveVal()
  restore_df_sequence_causality <- reactiveVal()

  # Input data upload
  df_upload <- reactive({
    data <- NULL
    if (input$data_input == 1) {
      if (is.null(input$upload)) {
        # From uploaded file
        data <- NULL
      } else {
        isolate({
          print(input$upload$type)
          print(input$upload$datapath)
          tryCatch({
          if ((tolower(input$upload$type) == "text/csv") ||
              (tolower(getExtension(input$upload$datapath)) == "csv")) {
            df_input_list <- lapply(input$upload$datapath, read.csv)
          } else if ((tolower(input$upload$type) == "text/tsv") ||
                (tolower(getExtension(input$upload$datapath)) == "tsv")) {
            df_input_list <- lapply(input$upload$datapath, read.tsv)
          } else if (grepl("xls", tolower(input$upload$type)) ||
                     grepl("excel", tolower(input$upload$type))) {
            df_input_list <- lapply(input$upload$datapath, read_excel)
          }
            names(df_input_list) <- gsub(input$upload$name, pattern="\\..*", replacement="")

            df_input <- bind_rows(df_input_list, .id = "group.file.chromo")

            data <- df_input
          }, error=function(cond){
            showNotification(paste("Could not open file because: ", cond), type = "error")
          })
        })
      }
    } else if (input$data_input == 2) {
      # From URL
      if (input$URL == "") {
        data <- NULL
      } else if (url.exists(input$URL) == FALSE) {
        data <- NULL
      } else {
        data <- read.csv(input$URL)
      }
    } else if (input$data_input == 3) {
      data <- example.data
    } else if (input$data_input == 4) {
      data <- example.data.huge
    }
    segmentation_raw(NULL)
    return(data)
  })

  # Display Data that has been uploaded
  output$data_uploaded <- renderDataTable({
    df_normalized()
  })

  output$fileUploaded <- reactive({
    return(!is.null(df_normalized()))
  })

  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)

  # Normalization
  df_normalized <- reactive ({
    cols_normalize <- input$cols_normalize
    if (input$norm_type != "none" && length(cols_normalize) > 0) {
      if (!all(sapply(df_selected()[,cols_normalize], is.numeric))) {
        showNotification("Only numeric columns can be normalized.", type = "error")
        norm_temp <- df_selected()
      } else if (input$norm_type == "max") {
        norm_temp <-
          df_selected() %>% group_by(!!sym(input$particle_vars)) %>%
          mutate_at(cols_normalize, normalize.max)
      } else if (input$norm_type == "min") {
        norm_temp <-
          df_selected() %>% group_by(!!sym(input$particle_vars)) %>%
          mutate_at(cols_normalize, normalize.min)
      } else if (input$norm_type == "median") {
        norm_temp <-
          df_selected() %>% group_by(!!sym(input$particle_vars)) %>%
          mutate_at(cols_normalize, normalize.median)
      } else if (input$norm_type == "zero_one") {
        norm_temp <-
          df_selected() %>% group_by(!!sym(input$particle_vars)) %>%
          mutate_at(cols_normalize, normalize.01)
      } else if (input$norm_type == "z-score") {
        norm_temp <-
          df_selected() %>% group_by(!!sym(input$particle_vars)) %>%
          mutate_at(cols_normalize, normalize.z)
      }
    } else {
      norm_temp <- df_selected()
    }

    if (is.null(norm_temp)) {
      return(NULL)
    }
    # Add velocities if selected
    coords <- input$coord_vars
    if (length(coords) != 2) {
      showNotification("Only tested with 2D analysis.", type = "warning")
    }
    if (input$calc_veloc) {
      norm_temp <- norm_temp %>% group_by(!!sym(input$particle_vars), !!sym(input$grouping_vars)) %>%
                  do(mutate(., cal.speed =
                          calculate.velocity(.,
                                        coords = coords),
                         ang.speed =
                          calc.angular.speed(.,
                                        coords = coords)))
    }
    if (input$calc_more) {
      norm_temp <- norm_temp %>% group_by(!!sym(input$particle_vars), !!sym(input$grouping_vars)) %>%
                  do(mutate(., cum.dist =
                          calculate.cumdist(.,
                                        coords = coords)))
    }
    norm_temp <- norm_temp %>%
                        na.omit() %>% ungroup()
    norm_temp <- as.data.frame(norm_temp)
    return(norm_temp)
  })

  df_process <- reactive ({
    df <- df_normalized()

    if (input$particle_vars == "none") {
      df <- df %>% mutate(particle = factor(1))
    } else {
      names(df)[names(df) == input$particle_vars] <- "particle"
    }
    if (input$grouping_vars == "none") {
      df <- df %>% mutate(group = factor(1))
    } else {
      names(df)[names(df) == input$grouping_vars] <- "group"
    }

    names(df)[names(df) == input$time_vars] <- "frame"

    return(df)
  })

  df_selected <- reactive({
    df <- df_upload()

    columns_select <- c(input$time_vars,
                        input$coord_vars,
                        input$particle_vars,
                        input$grouping_vars)

    if (!is.empty(columns_select) && input$show_selected) {
      df <- df[, columns_select]
    }

    if (input$rotate_coords_seg) {
      library(spdep)
      df <-
        as.data.frame(rotate.axes.maximum(df, input$coord_vars, input$particle_vars))
      detach("package:spdep", unload=TRUE)
    }

    if (input$exclude_cells && !is.empty(input$cells_to_remove)) {
      df <- df %>% filter(!(!!sym(input$particle_vars) %in% input$cells_to_remove))
    }

    if (input$norm_time == "last") {
      df <- df %>% group_by(!!sym(input$particle_vars)) %>%
        mutate(!!sym(input$time_vars) := !!sym(input$time_vars) - max(!!sym(input$time_vars))) %>%
        ungroup()
      df <- as.data.frame(df)
    }

    if (input$norm_time == "last_norm") {
      df <- df %>% group_by(!!sym(input$particle_vars)) %>%
        mutate(!!sym(input$time_vars) := !!sym(input$time_vars) - max(!!sym(input$time_vars))) %>%
        ungroup() %>% mutate(!!sym(input$time_vars) := !!sym(input$time_vars) - min(!!sym(input$time_vars)))
      df <- as.data.frame(df)
    }

    return(df)
  })

  observe({
    if (is.null(df_upload())) {
      return()
    }
    var_names  <- names(df_upload())
    var_list <- c("none", var_names)

    input_names_cols <- c(
      "cols_normalize",
      "time_vars",
      "coord_vars",
      "particle_vars",
      "grouping_vars"
    )

    for (u in input_names_cols) {
      updateSelectizeInput(session, u,
                           choices = var_list)
    }
  })
  
  observe({
    if (input$adjustcolors == "1") {
      plotColor <<- styles_color_Default
    } else if (input$adjustcolors == "2") {
      plotColor <<- gsub("\\s","", strsplit(input$user_color_list,",")[[1]])
    }
  })

  observe({
    if (is.null(df_upload())) {
      return()
    }
    var_names  <- names(df_upload())
    var_list <- c("none", var_names)

    input_names_cols_calc <- c(
      "cols_individual_plot",
      "df_vars_motifs",
      "segment_variables",
      "density_variables",
      "df_query_causality"
    )

    if (input$calc_veloc) {
      var_list <- c(var_list,
                    "ang.speed",
                    "cal.speed")
    }

    if (input$calc_more) {
      var_list <- c(var_list,
                    "cum.dist")
    }

    for (u in input_names_cols_calc) {
      updateSelectizeInput(session, u,
                           choices = var_list)
    }

    # Update values when bookmarked state is available
    if(!is.empty.or.except(restore_segment_variables())) {
      updateSelectizeInput(session, "segment_variables",
                           selected = restore_segment_variables())
    }

    if(!is.empty.or.except(restore_df_vars_motifs())) {
      updateSelectizeInput(session, "df_vars_motifs",
                           selected = restore_df_vars_motifs())
    }

    if(!is.empty.or.except(restore_df_query_causality())) {
      updateSelectizeInput(session, "df_query_causality",
                           selected = restore_df_query_causality())
    }

    # Find the coordinate variables
    coords <-
      intersect(var_list,
                c("time", "Time", "frame", "Frame", "slice", "Slice"))

    updateSelectizeInput(session, "time_vars",
                         selected = coords)

    # Find the coordinate variables
    coords <-
      intersect(
        var_list,
        c(
          "x",
          "y",
          "z",
          "coord.x",
          "coord.y",
          "coord.z",
          "X",
          "Y",
          "Z",
          "COORD.X",
          "COORD.Y",
          "COORD.Z"
        )
      )

    updateSelectizeInput(session, "coord_vars",
                         selected = coords)

    # Find the label variable
    label <-
      intersect(var_list,
                c("label", "group", "strain", "Label", "Group", "Strain", "group.file.chromo"))

    updateSelectizeInput(session, "grouping_vars",
                         selected = label)

    # Find the particle variable
    particle <-
      intersect(var_list, c("particle", "cell", "Particle", "Cell"))

    updateSelectizeInput(session, "particle_vars",
                         selected = particle)
  })

  observe({
    if (input$particle_vars != "none" && !is.null(df_upload())) {
      particles <-
        df_normalized() %>% select(!!sym(input$particle_vars)) %>% distinct()
      input_names_part <- c("part_individual_plot")
      for (u in input_names_part) {
        updateSelectizeInput(session, u,
                             choices = particles[, 1])
      }
    }
  })

  observe({
    if (input$particle_vars != "none" && !is.null(df_upload())) {
      particles <-
        df_upload() %>% select(!!sym(input$particle_vars)) %>% distinct()
      input_names_part <- c("cells_to_remove")
      for (u in input_names_part) {
        updateSelectizeInput(session, u,
                             choices = particles[, 1])
      }
    }
  })

  observe({
    if (input$grouping_vars != "none" && !is.null(df_upload())) {
      groups <-
        df_normalized() %>% select(!!sym(input$grouping_vars)) %>% distinct()
      input_names_group <- c("df_query_motifs")
      for (u in input_names_group) {
        updateSelectizeInput(session,
                             u,
                             choices = c("none", groups[, 1]),
                             selected = "none")
      }

      if(!is.empty.or.except(restore_df_sequence_motifs())) {
        selected_motifs <- restore_df_sequence_motifs()
      } else {
        selected_motifs <- groups[1, 1]
      }

      if(!is.empty.or.except(restore_df_sequence_causality())) {
        selected_causality <- restore_df_sequence_causality()
      } else {
        selected_causality <- groups[1, 1]
      }
      updateSelectizeInput(session,
                           "df_sequence_motifs",
                           choices = groups[, 1],
                           selected = selected_motifs)
      updateSelectizeInput(session,
                           "df_sequence_causality",
                           choices = groups[, 1],
                           selected = selected_causality)
    }
  })


  segmentation_raw <- reactiveVal()
  observeEvent(input$run_segmentation, {
    validate(
      need(input$segment_variables, "\n\n\n\nPlease, select at least one column to segment"
      )
    )
    # Parsing reactive values
    df_in <- df_process()
    seg_method <- input$seg_method
    segment_variables <- input$segment_variables
    sampl_freq <- input$sampl_freq

    # Creating process handler
    progress <- Progress$new(session)
    progress$set(message = "Computing segmentation", value = 0)

    future(seed=TRUE, {
      updateProgress <- function(detail = NULL) {
        progress$inc(amount = 1/4, detail = detail)
      }

      if(seg_method == "lavielle") {
        segments <- segment.timeseries(df_in,
                    coords = segment_variables,
                    sampl_freq,
                    updateProgress = NULL)
      } else {
        segments <- segment.timeseries.multiv(df_in,
                    coords = segment_variables,
                    sampling.time = sampl_freq,
                    updateProgress = NULL)
      }


      progress$close()
      segments
    }) %...>% segmentation_raw()
  })


  segmentation_calc <- reactiveVal()
  observe({
    segments <- clusterize.segments(req(segmentation_raw()[[1]]),
                                    input$grouping_vars,
                                    input$particle_vars,
                                    segments = input$seg_range)
    # segments <- reapply.clustering(segments)
    clust <- create.df.clusters(segments, segmentation_raw()[[2]])
    segmentation_clusters(clust)
    segmentation_calc(segments)
  })

  plot_individual <- reactive({
    if (input$tabsegment == "Individual" &&
        !is.empty(input$cols_individual_plot) &&
        !is.empty(input$part_individual_plot)) {

      plt_df <- df_normalized()
      plt_df$part.chromo <- as.factor(plt_df[[input$particle_vars]])
      if (!is.empty(input$part_individual_plot)) {
        plt_df <-
          plt_df %>% filter(part.chromo == input$part_individual_plot)
      }
      p <- plt_df %>%
        pivot_longer(
          cols = input$cols_individual_plot,
          names_to = "variable",
          values_to = "values"
        ) %>%
        ggplot(aes(
          !!sym(input$time_vars),
          values,
          color = !!sym(input$particle_vars)
        )) +
        facet_grid(variable * part.chromo ~  ., scales = "free") +
        geom_path(aes(group = !!sym(input$particle_vars))) +
        theme_bw() +
        scale_color_manual(values=plotColor) +
        scale_fill_manual(values=plotColor) +
        labs(x = "Time") +
        theme(legend.position = "none")

      return(p)
    }
  })

  plot_segments <- reactive({
    validate(
      need(segmentation_calc(), "\n\n\n\nPlease, run 'Calculate Segmentation' to continue."
      )
    )
    p <- segmentation.plot(segmentation_calc())
    return(p)
  })

  plot_segments_ly <- reactive({
    validate(
      need(segmentation_calc(), "\n\n\n\nPlease, run 'Calculate Segmentation' to continue."
      )
    )
    p <- segmentation.plotly(segmentation_calc())
    return(p)
  })

  plot_spectrum <- reactive({
    validate(
      need(spectrum_global(), "")
    )

    val <- spectrum_global()["spec.f"]
    if (!input$freq_spectrogram) {
      val <- (1 / val) / 60
    }

    updateSliderInput(session, "spec_range", value = c(min(val)*0.1, max(val)*0.6),
          min = min(val), max = max(val))

    p <- spectral.plot(
      spectrum_global(),
      input$sampl_freq,
      input$freq_spectrogram,
      input$spectrum_individual,
      input$spec_ticks
    )
    return(p) # in Plotly format
  })

  plot_spectrum_segments <- reactive({
    validate(
      need(spectrum_segments(), "\nPlease, run 'Calculate Spectrum' to continue."
      )
    )
    p <- spectral.plot.segments(
      spectrum_segments(),
      input$sampl_freq,
      input$freq_spectrogram
    )
    return(p) # in Plotly format
  })

  plot_heatmap <- reactive({
    validate(
      need(spectrogram_global(), "")
    )
    p <- heatmap.plot(
      spectrogram_global(),
      input$freq_spectrogram
    )
    return(p) # in Plotly format
  })

  plot_heatmap_segment <- reactive({
    validate(
      need(spectrogram_segments(), "\nPlease, run 'Calculate Spectrum' to continue."
      )
    )
    p <- heatmap.segmented.plot(
      spectrogram_segments(),
      input$freq_spectrogram
    )
    return(p) # in Plotly format
  })

  plot_velocities_heat <- reactive({
    validate(
      need(velocities_group(), "")
    )
    p <- heatmap.veloc.plot(
      velocities_group(),
      input$density_variables
    )
    return(p) # in Plotly format
  })

  plot_velocities_heat_segment <- reactive({
    validate(
      need(velocities_segment(), "\nPlease, run 'Calculate Densities' to continue."
      )
    )
    p <- heatmap.veloc.segmented.plot(
      velocities_segment(),
      input$density_variables
    )
    return(p) # in Plotly format
  })

  segment_significant <- reactive({
    t <-
      segmentation.significance(segmentation_calc(),
                                input$signif_segment_test,
                                input$signif_segment_adj)
    return(t)
  })

  adj_matrix_show <- reactive({
    validate(
      need(causality_discovered(), "Please, run 'Causality Analysis' to continue."
      )
    )
    return(causality_discovered())
  })

  segment_clustering <- reactive({
    t <- segmentation.summary(segmentation_calc())
    return(t)
  })

  segment_clustering_grouped <- reactive({
    t <- segmentation.grouped.summary(segmentation_calc())
    return(t)
  })

  motifs_discovery <- reactiveVal()
  observe({
    if (input$df_vars_motifs == "none") {
      return(NULL)
    }

    if (!all(sapply(df_normalized()[,input$df_vars_motifs], is.numeric))) {
        showNotification("Only numeric columns can be analyzed.", type = "error")
        return(NULL)
    }

    df <- df_normalized()
    seqs <-
      df[df[[input$grouping_vars]] == input$df_sequence_motifs, input$df_vars_motifs]

    time <-
      df[df[[input$grouping_vars]] == input$df_sequence_motifs, input$time_vars]
    groups <-
      df[df[[input$grouping_vars]] == input$df_sequence_motifs, input$grouping_vars]
    particles <-
      df[df[[input$grouping_vars]] == input$df_sequence_motifs, input$particle_vars]

    motifs <-
      motif.discovery(
        seqs,
        win = input$motif_window,
        pct = input$motif_sample_pct,
        thr = input$motif_correlation,
        time = time,
        group = groups,
        particles = particles,
        nmotifs = input$motif_amount
      )
    motifs_discovery(motifs)
  })

  motifs_discovery_segmented <- reactiveVal()

  observeEvent(input$run_motifs, {
    if(is.empty(segmentation_clusters())) {
      showNotification("You have to run a segmentation, first!", type = "error")
    }
    req(segmentation_clusters())
    if (input$df_vars_motifs == "none") {
      return(NULL)
    }

    if (!all(sapply(df_normalized()[,input$df_vars_motifs], is.numeric))) {
        showNotification("Only numeric columns can be analyzed.", type = "error")
        return(NULL)
    }

    progress <- Progress$new(session)
    progress$set(message = "Computing Motifs", value = 0)

    segments.all <- segmentation_clusters()

    # Temporally retrieving reactive values for future
    df_vars_motifs <- input$df_vars_motifs
    df_sequence_motifs <- input$df_sequence_motifs
    motif_window <- input$motif_window
    motif_sample_pct <- input$motif_sample_pct
    motif_correlation <- input$motif_correlation
    motif_amount <- input$motif_amount

    future(seed=TRUE, {
      motif.list <- motifs.per.segment(segments.all,
                                     df_vars_motifs,
                                     df_sequence_motifs,
                                     NULL,
                                     motif_window,
                                     motif_sample_pct,
                                     motif_correlation,
                                     motif_amount)
      progress$close()
      motif.list
    }) %...>% motifs_discovery_segmented()
  })

  causality_discovered <- reactiveVal()

  cor_group <- reactiveVal()

  observeEvent(input$run_causal, {
    req(input$df_query_causality)
    req(input$df_sequence_causality)
    req(df_process())

    if (!all(sapply(df_process()[,input$df_query_causality], is.numeric))) {
        showNotification("Only numeric columns can be analyzed.", type = "error")
        return()
    }

    segments.all <- segmentation_clusters()
    df <- df_process()

    progress <- Progress$new(session)
    progress$set(message = "Computing Causality", value = 0)

    # Temporally retrieving reactive values for future
    df_query_causality <- input$df_query_causality
    pval_corr_range <- input$pval_corr_range
    lags_corr_range <- input$lags_corr_range
    signif_causal_adj <- input$signif_causal_adj
    df_sequence_causality <- input$df_sequence_causality

    cor_group(calc.cor(df,
                          df_query_causality,
                          df_sequence_causality))

    causality.list <- causality.global(df,
                          df_query_causality,
                          pval_corr_range,
                          lags_corr_range,
                          signif_causal_adj,
                          df_sequence_causality,
                          progress = progress)

    progress$close()

    causality_discovered(causality.list)
  })

  spectrum_global <- reactiveVal()
  spectrogram_global <- reactiveVal()
  spectrum_segments <- reactiveVal()
  spectrogram_segments <- reactiveVal()

  observeEvent(input$run_spectrum, {
    progress <- Progress$new(session)
    progress$set(message = "Computing Spectrum", value = 0)

    df <- df_process()
    segments.all <- segmentation_clusters()

    # Temporally retrieving reactive values for future
    coord_vars <- input$coord_vars
    sampl_freq <- input$sampl_freq

    future(seed=TRUE, {

      updateProgress <- function(detail = NULL) {
        progress$inc(amount = 1/5, detail = detail)
      }

      df_freqs <- spectrum.per.segment(
                  segments.all,
                  coord_vars[1],
                  sampl_freq,
                  updateProgress = updateProgress
                )
      progress$close()
      df_freqs

    }) %...>% spectrum_segments()

    progress.2 <- Progress$new(session)
    progress.2$set(message = "Computing Spectrograms", value = 0)

    if(is.empty(segmentation_clusters())) {
      showNotification("You have to run a segmentation, first!", type = "error")
    }
    req(segmentation_clusters())

    future(seed=TRUE, {

      updateProgress <- function(detail = NULL) {
        progress.2$inc(amount = 1/5, detail = detail)
      }

      df_freqs <- spectrogram.per.segment(
                    segments.all,
                    coord_vars,
                    updateProgress = updateProgress
                  )
      progress.2$close()
      df_freqs

    }) %...>% spectrogram_segments()

    progress.3 <- Progress$new(session)
    progress.3$set(message = "Computing Global Spectrum", value = 0)

    future(seed=TRUE, {

      updateProgress <- function(detail = NULL) {
        progress.3$inc(amount = 1/5, detail = detail)
      }

      df_freqs <- spectrum.global(
                    df,
                    coord_vars,
                    sampl_freq,
                    updateProgress = updateProgress
                  )
      progress.3$close()
      df_freqs

    }) %...>% spectrum_global()

    progress.4 <- Progress$new(session)
    progress.4$set(message = "Computing Global Spectrogram", value = 0)

    future(seed=TRUE, {

      updateProgress <- function(detail = NULL) {
        progress.4$inc(amount = 1/5, detail = detail)
      }
      df_freqs <- spectrogram.global(
                    segments.all,
                    coord_vars,
                    updateProgress = updateProgress
                  )
      progress.4$close()
      df_freqs

    }) %...>% spectrogram_global()
  })

  plot_motifs <- reactive({
    return(req(motifs_discovery()))
  })

  plot_motifs_segmented <- reactive({
    return(req(motifs_discovery_segmented()))
  })

  spectrum_significance <- reactive({
    validate(
      need(spectrum_global(), "Please, run 'Calculate Spectrum' to continue."
      )
    )
    spec_signif <- spectrum.significance(
      spectrum_global(),
      input$spec_range
    )
    return(spec_signif)
  })

  spectrum_significance_segment <- reactive({
    validate(
      need(spectrum_global(), "Please, run 'Calculate Spectrum' to continue."
      )
    )
    spec_signif <- spectrum.significance.segment(
      spectrum_segments(),
      input$spec_range
    )
    return(spec_signif)
  })

  plot_morlet <- reactive({
    plt_ml <- df_normalized()
    plt_ml <-
      plt_ml %>% filter(!!sym(input$particle_vars) == input$part_individual_plot)
    plt <-
      morlet.plot(plt_ml, input$cols_individual_plot[1],
                  sig = input$ind_spec_level,
                  ylab = input$cols_individual_plot)
    return(plt)
  })

  velocities_group <- reactiveVal()
  velocities_segment <- reactiveVal()
  segmentation_clusters <- reactiveVal()

  observeEvent(input$run_densities, {
    if (input$density_variables == "none") {
      return(NULL)
    }

    df <- df_process()
    df.clu <- segmentation_clusters()

    # Temporally retrieving reactive values for future
    particle_vars <- input$particle_vars
    density_variables <- input$density_variables
    vel_ma <- input$vel_ma

    future(seed=TRUE, {
      df_vt.smooth <- calc.densities(df,
                                    density_variables,
                                    vel_ma)
      df_vt.smooth

    }) %...>% velocities_group()

    progress <- Progress$new(session)
    progress$set(message = "Computing Velocities", value = 0)

    if(is.empty(segmentation_clusters())) {
      showNotification("You have to run a segmentation, first!", type = "error")
    }
    req(segmentation_clusters())

    future(seed=TRUE, {

      updateProgress <- function(detail = NULL) {
        progress$inc(amount = 1/2, detail = detail)
      }

      df_seg_vels <- density.per.segment(df.clu,
                      density_variables,
                      vel_ma,
                      updateProgress = updateProgress)
      progress$close()
      df_seg_vels

    }) %...>% velocities_segment()
  })

  msd_group <- reactiveVal()

  msd_segment <- reactiveVal()

  observeEvent(input$run_msd, {
    segments.all <- segmentation_clusters()
    df <- df_process()

    # Temporally retrieving reactive values for future
    coord_vars <- input$coord_vars
    displ_k <- input$displ_k

    msd_group(calc.msd(df,
                          coord_vars,
                          displ_k))

    if(is.empty(segmentation_clusters())) {
      showNotification("You have to run a segmentation, first!", type = "error")
    }
    req(segmentation_clusters())

    msd_segment(calc.msd.segment(segments.all,
                                  coord_vars,
                                  displ_k))
  })

  plot_velocities <- reactive({
    validate(
      need(velocities_group(), "")
    )
    p <- plot.velocities(
      velocities_group(),
      input$densities_individual,
      input$density_variables
    )
    return(p)
  })

  plot_velocities_segment <- reactive({
    validate(
      need(velocities_segment(), "\nPlease, run 'Calculate Densities' to continue."
      )
    )
    p <- plot.velocities.segment(
      velocities_segment(),
      input$density_variables
    )
    return(p)
  })

  plot_msd <- reactive({
    validate(
      need(msd_group(), "")
    )
    p <- plot.msd(
      msd_group()
    )
    return(p)
  })

  plot_msd_segment <- reactive({
    validate(
      need(msd_segment(), "\nPlease, run 'Calculate MSD' to continue."
      )
    )
    p <- plot.msd.segment(
      msd_segment()
    )
    return(p)
  })

  plot_cor <- reactive({
    validate(
      need(cor_group(), "\nPlease, run 'Causality Analysis' to continue."
      )
    )
    return(chart.Correlation(cor_group(), histogram=TRUE, pch="+"))
  })

  plot_causality_pc <- reactive({
    validate(
      need(causality_discovered(), "\nPlease, run 'Causality Analysis' to continue."
      )
    )
    causality.plot(causality_discovered()[[1]],
                      input$perc_corr_range,
                      input$df_query_causality)
  })

  plot_causality_vlte <- reactive({
    validate(
      need(causality_discovered(), "\nPlease, run 'Causality Analysis' to continue."
      )
    )
    causality.plot(causality_discovered()[[2]],
                      input$perc_corr_range,
                      input$df_query_causality)
  })

  plot_causality_vlgc <- reactive({
    validate(
      need(causality_discovered(), "\nPlease, run 'Causality Analysis' to continue."
      )
    )
    causality.plot(causality_discovered()[[3]],
                      input$perc_corr_range,
                      input$df_query_causality)
  })

  velocities_significance <- reactive({
    validate(
      need(velocities_group(), "Please, run 'Calculate Densities' to continue."
      )
    )
    significant.velocities(velocities_group())
  })

  velocities_significance_segment <- reactive({
    validate(
      need(velocities_group(), "Please, run 'Calculate Densities' to continue."
      )
    )
    significant.velocities.segment(velocities_segment())
  })

  plot_width <- reactive ({
    input$plot_width
  })
  plot_height <- reactive ({
    input$plot_height
  })
  plot_height_spec <- reactive ({
    input$plot_height * 3
  })
  plot_height_motif <- reactive ({
    input$plot_height * 1.75
  })

  output$plot_segments <-
    renderPlotly({
      plot_segments_ly()
    })

  output$plot_spectrum <-
    renderPlotly({
      ggplotly(plot_spectrum(), width = 600)
    })

  output$plot_spectrum_segments <-
    renderPlotly({
      ggplotly(plot_spectrum_segments(), width = 600)
    })

  output$plot_heatmap <-
    renderPlot(width = plot_width, height = plot_height, {
      plot_heatmap()
    })

  output$plot_heatmap_segment <-
    renderPlot(width = plot_width, height = plot_height, {
      plot_heatmap_segment()
    })

  output$plot_individual <-
    renderPlotly({
      ggplotly(plot_individual(), width = 600)
    })

  output$table_segments <- renderDataTable({
    segment_clustering()
  })

  output$table_segments_grouped <- renderDataTable({
    segment_clustering_grouped()
  })

  output$table_segments_significant <- renderPrint({
    segment_significant()
  })

  output$adj_matrix <- renderPrint({
    adj_matrix_show()
  })

  output$plot_motifs_a <-
    renderPlot(width = plot_width, height = plot_height, {
      plot(tsmp::as.matrixprofile(plot_motifs()[[1]]))
    })

  output$plot_motifs_b <-
    renderPlot(width = plot_width, height = plot_height, {
      plot(tsmp::motifs(plot_motifs()[[1]]))
    })

  output$plot_motifs_c <-
    renderPlot(width = plot_width, height = plot_height, {
      plot(tsmp::discords(plot_motifs()[[1]]))
    })

  output$plot_morlet <-
    renderPlot(width = plot_width, height = plot_height, {
      plot_morlet()
    })

  output$plot_velocities <-
    renderPlotly({
      ggplotly(plot_velocities(), width = 600)
    })

  output$plot_velocities_segment <-
    renderPlotly({
      ggplotly(plot_velocities_segment(), width = 600)
    })

  output$plot_velocities_heat <-
    renderPlot(width = plot_width, height = plot_height, {
      plot_velocities_heat()
    })

  output$plot_velocities_heat_segment <-
    renderPlot(width = plot_width, height = plot_height, {
      plot_velocities_heat_segment()
    })

  output$plot_msd <-
    renderPlotly({
      ggplotly(plot_msd(), width = 600)
    })

  output$plot_msd_segment <-
    renderPlotly({
      ggplotly(plot_msd_segment(), width = 600)
    })

  output$plot_cor <-
    renderPlot(width = plot_width, height = plot_height, {
      plot_cor()
    })

  output$summ_spectrum <- renderPrint({
    spectrum_significance()
  })

  output$summ_spectrum_segment <- renderPrint({
    spectrum_significance_segment()
  })

  output$summ_velocities <- renderPrint({
    velocities_significance()
  })

  output$summ_velocities_segment <- renderPrint({
    velocities_significance_segment()
  })

  output$tabs_motifs <- renderUI({
    seg.toplot <- req(plot_motifs_segmented())
    range <- c(1, 180)
    for (i in 1:length(seg.toplot)) {
      local({
        my_i <- i

        tab_s_motif <-
          paste("download_segment_motifs_", my_i, sep = "")

        output[[tab_s_motif]] <- downloadHandler(
          filename <- function() {
            paste("ChroMo_", Sys.time(), "_Segment_", i, "_Motifs.json", sep = "")
          },
          content <- function(file) {
            motifs_tab_download <- seg.toplot[[my_i]]
            write(rjson::toJSON(motifs_tab_download, indent=1), file)
          }
        )

        tab_s_motif <-
          paste("download_segment_motifs_csv", my_i, sep = "")

        output[[tab_s_motif]] <- downloadHandler(
          filename <- function() {
            paste("ChroMo_", Sys.time(), "_Segment_", i, "_Motifs.csv", sep = "")
          },
          content <- function(file) {
            motifs_tab_download <- seg.toplot[[my_i]][[2]]
            write.csv(motifs_tab_download, file, row.names = FALSE)
          }
        )

        plt_discord <-
          paste("plot_segment_discord_", my_i, sep = "")

        output[[plt_discord]] <-
          renderPlot(width = plot_width, height = plot_height, {
            plot(tsmp::discords(seg.toplot[[my_i]][[1]]))
          })

        plt_dur_discord <-
          paste("plot_location_discord_", my_i, sep = "")

        output[[plt_dur_discord]] <-
          renderPlotly({
            ggplotly(locations.motifs.plot(seg.toplot[[my_i]][[3]]), width = 600, height = 200)
          })

        plt_motif <- paste("plot_segment_motifs_", my_i, sep = "")

        output[[plt_motif]] <- renderPlot(width = plot_width, height = plot_height, {
          plot(tsmp::motifs(seg.toplot[[my_i]][[1]]))
        })

        plt_dur_motif <-
          paste("plot_location_motifs_", my_i, sep = "")

        output[[plt_dur_motif]] <-
          renderPlotly({
            ggplotly(locations.motifs.plot(seg.toplot[[my_i]][[2]]), width = 600, height = 200)
          })
      })
    }
    plot_output_list <- lapply(1:length(seg.toplot), function(i) {
      tabPanel(paste("Segment ", i, sep = ""),
               hr(),
               downloadButton(paste("download_segment_motifs_", i, sep = ""), "Download JSON"),
               downloadButton(paste("download_segment_motifs_csv", i, sep = ""), "Download CSV"),
               plotOutput(paste("plot_segment_motifs_", i, sep = "")),
               plotlyOutput(paste("plot_location_motifs_", i, sep = ""), height = "auto"),
               plotOutput(paste("plot_segment_discord_", i, sep = "")),
               plotlyOutput(paste("plot_location_discord_", i, sep = "")), height = "auto")
    })

    do.call(tabsetPanel, plot_output_list)
  })

  output$plot_vlgc_corr <-
    renderForceNetwork(plot_causality_vlgc())

  output$plot_pc_corr <-
    renderForceNetwork(plot_causality_pc())

  output$plot_vlte_corr <-
    renderForceNetwork(plot_causality_vlte())

  output$downSegm <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_Segments.pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width / 72, height = input$plot_height / 72)
      plot(plot_segments())
      dev.off()
    },
    contentType = "application/pdf"
  )

  output$downSpec <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_Spectrum.pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width / 72, height = input$plot_height / 72)
      plot(plot_spectrum())
      dev.off()
    },
    contentType = "application/pdf"
  )

  output$downSpecTable <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_Spectrum.csv", sep = "")
    },
    content <- function(file) {
      spectrum_tab_download <- spectrum_global()[c("particle",
                                                 "group",
                                                 "spec.s",
                                                 "spec.f")]

      write.csv(spectrum_tab_download, file, row.names = FALSE)
    }
  )

  output$downSpecSeg <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_Spectrum.pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width / 72, height = input$plot_height / 72)
      plot(plot_spectrum_segments())
      dev.off()
    },
    contentType = "application/pdf"
  )

  output$downSpecH <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_Spectrum.pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width / 72, height = input$plot_height / 72)
      plot(plot_heatmap())
      dev.off()
    },
    contentType = "application/pdf"
  )

  output$downSpecSegTable <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_Spectrum_Segments.csv", sep = "")
    },
    content <- function(file) {
      spectrum_tab_download <- spectrum_segments()[c("particle",
                                                     "group",
                                                     "segmentchromo",
                                                     "frame.init",
                                                     "frame.seg",
                                                     "cluster",
                                                     "spec.s",
                                                     "spec.f")]

      write.csv(spectrum_tab_download, file, row.names = FALSE)
    }
  )

  output$downSpecSegH <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_Spectrum.pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width / 72, height = input$plot_height / 72)
      plot(plot_heatmap_segment())
      dev.off()
    },
    contentType = "application/pdf"
  )

  output$downVel <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_Velocities.pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width / 72, height = input$plot_height / 72)
      plot(plot_velocities())
      dev.off()
    },
    contentType = "application/pdf"
  )

  output$downVelTable <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_", input$density_variables, ".csv", sep = "")
    },
    content <- function(file) {
      density.tab.dw <- velocities_group()[c("particle",
                                             "group",
                                             "vel.ma")]
      names(density.tab.dw)[names(density.tab.dw) == 'vel.ma'] <- input$density_variables
      write.csv(density.tab.dw, file, row.names = FALSE)
    }
  )

  output$downVelSeg <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_Velocities.pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width / 72, height = input$plot_height / 72)
      plot(plot_velocities_segment())
      dev.off()
    },
    contentType = "application/pdf"
  )

  output$downVelSegTable <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_", input$density_variables, "_Segments.csv", sep = "")
    },
    content <- function(file) {
      density.tab.dw <- velocities_segment()[c("particle",
                                             "group",
                                             "segmentchromo",
                                             "vel.ma")]
      names(density.tab.dw)[names(density.tab.dw) == 'vel.ma'] <- input$density_variables
      write.csv(density.tab.dw, file, row.names = FALSE)
    }
  )

  output$downVelH <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_Velocities.pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width / 72, height = input$plot_height / 72)
      plot(plot_velocities_heat())
      dev.off()
    },
    contentType = "application/pdf"
  )

  output$downVelSegH <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_Velocities.pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width / 72, height = input$plot_height / 72)
      plot(plot_velocities_heat_segment())
      dev.off()
    },
    contentType = "application/pdf"
  )

  output$downMSD <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_MSD.pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width / 72, height = input$plot_height / 72)
      plot(plot_msd())
      dev.off()
    },
    contentType = "application/pdf"
  )

  output$downMSDSeg <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_MSD.pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width / 72, height = input$plot_height / 72)
      plot(plot_msd_segment())
      dev.off()
    },
    contentType = "application/pdf"
  )

  output$downIndiv <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_MSD.pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width / 72, height = input$plot_height / 72)
      plot(plot_individual())
      dev.off()
    },
    contentType = "application/pdf"
  )

  output$downMP <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_MatrixProfile.pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width / 72, height = input$plot_height / 72)
      plot(as.matrixprofile(plot_motifs()))
      dev.off()
    },
    contentType = "application/pdf"
  )

  output$downMotif <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_Motifs.pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width / 72, height = input$plot_height / 72)
      plot(tsmp::motifs(plot_motifs()))
      dev.off()
    },
    contentType = "application/pdf"
  )

  output$downMotifTab <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_Motifs.json", sep = "")
    },
    content <- function(file) {
      motifs_tab_download <- motifs_discovery()
      write(rjson::toJSON(motifs_tab_download, indent=1), file)
    }
  )

  output$downCSVMotifTab <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_Motifs.csv", sep = "")
    },
    content <- function(file) {
      motifs_tab_download <- motifs_discovery()[[2]]
      write.csv(motifs_tab_download, file, row.names = FALSE)
    }
  )

  output$downDiscord <- downloadHandler(
    filename <- function() {
      paste("ChroMo_", Sys.time(), "_Discords.pdf", sep = "")
    },
    content <- function(file) {
      pdf(file, width = input$plot_width / 72, height = input$plot_height / 72)
      plot(tsmp::discords(plot_motifs()))
      dev.off()
    },
    contentType = "application/pdf"
  )

  output$downSpecgram <- downloadHandler(
    filename <- function() {data:image
      paste("ChroMo_", Sys.time(), "_Specgram.pdf", sep = "")
    },
    content <- function(file) {
      plot_morlet()
      dev.copy(pdf, file = file, width = 4, height = 4)
      dev.off()
    },
    contentType = "application/pdf"
  )

  output$downPC <- downloadHandler(
    filename <- function() {
      paste("ChroMo_Causation_PC.html", sep = "")
    },
    content <- function(file) {
      saveNetwork(plot_causality_pc() %>%
        htmlwidgets::prependContent(
          htmltools::tags$h1("Peter-Clark algorithm directed graph")
        ),
      file, selfcontained = TRUE)
    },
    contentType = "text/html"
  )

  output$downVLTE <- downloadHandler(
    filename <- function() {
      paste("ChroMo_Causation_VLTE.html", sep = "")
    },
    content <- function(file) {
      saveNetwork(plot_causality_vlte() %>%
        htmlwidgets::prependContent(
          htmltools::tags$h1("Variable Lag Transfer Entropy directed graph")
        ),
      file, selfcontained = TRUE)
    },
    contentType = "text/html"
  )

  output$downVLGC <- downloadHandler(
    filename <- function() {
      paste("ChroMo_Causation_VLGC.html", sep = "")
    },
    content <- function(file) {
      saveNetwork(plot_causality_vlgc() %>%
        htmlwidgets::prependContent(
          htmltools::tags$h1("Variable Lag Granger Causality directed graph")
        ),
      file, selfcontained = TRUE)
    },
    contentType = "text/html"
  )

  output$downCSV <- downloadHandler(
    filename <- function() {
      paste("ChroMo_data.csv", sep = "")
    },
    content <- function(file) {
      write.csv(df_normalized(), file, row.names = FALSE)
    }
  )

  output$count <- renderText({
    vals$count
  })

  session$onSessionEnded(function() {
    isolate(vals$count <- vals$count - 1)
  })

  onBookmark(function(state) {
    state$values$segmentation_raw <- segmentation_raw()
    state$values$motifs_discovery_segmented <- motifs_discovery_segmented()
    state$values$motifs_discovery <- motifs_discovery()
    state$values$causality_discovered <- causality_discovered()
    state$values$cor_group <- cor_group()
    state$values$spectrum_global <- spectrum_global()
    state$values$spectrogram_global <- spectrogram_global()
    state$values$spectrum_segments <- spectrum_segments()
    state$values$spectrogram_segments <- spectrogram_segments()
    state$values$velocities_group <- velocities_group()
    state$values$velocities_segment <- velocities_segment()
    state$values$segmentation_clusters <- segmentation_clusters()
    state$values$segmentation_calc <- segmentation_calc()
    state$values$msd_group <- msd_group()
    state$values$msd_segment <- msd_segment()
  })

  # Read values from state$values when we restore
  onRestore(function(state) {
    segmentation_raw(state$values$segmentation_raw)
    motifs_discovery_segmented(state$values$motifs_discovery_segmented)
    motifs_discovery(state$values$motifs_discovery)
    causality_discovered(state$values$causality_discovered)
    cor_group(state$values$cor_group)
    spectrum_global(state$values$spectrum_global)
    spectrogram_global(state$values$spectrogram_global)
    spectrum_segments(state$values$spectrum_segments)
    spectrogram_segments(state$values$spectrogram_segments)
    velocities_group(state$values$velocities_group)
    velocities_segment(state$values$velocities_segment)
    segmentation_clusters(state$values$segmentation_clusters)
    segmentation_calc(state$values$segmentation_calc)
    msd_group(state$values$msd_group)
    msd_segment(state$values$msd_segment)
    restore_segment_variables(state$input$segment_variables)
    restore_df_vars_motifs(state$input$df_vars_motifs)
    restore_df_query_causality(state$input$df_query_causality)
    restore_df_sequence_motifs(state$input$df_sequence_motifs)
    restore_df_sequence_causality(state$input$df_sequence_causality)
  })

}