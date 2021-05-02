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

options(shiny.maxRequestSize = 5 * 1024 ^ 2) # maximum 5 mb
library(shiny)
library(shinycssloaders)
library(promises)
library(future)
library(ipc)
library(forecast)
library(dplyr)
library(ggplot2)
library(dplR)
library(lmtest)
library(purrr)
library(tidyr)
library(reshape2)
library(tseries)
library(segclust2d)
library(segmenTier)
library(mclust)
library(lattice)
library(TTR)
library(egg)
library(readxl)
library(RCurl)
library(shinythemes)
library(ggpubr)
library(htmlwidgets)
library(plotly)
library(rapportools)
library(tie)
# For causality analysis
library(pcalg)
library(VLTimeCausality)
library(network)
library(sna)
library(networkD3)
library(igraph)
library(GGally)
library(PerformanceAnalytics)

plan(sequential)

source("chromo_functions.R")
example.data <- read.csv("../examples/example.csv")
example.data.huge <- read.csv("../examples/example_huge.csv")
vals <- reactiveValues(count = 0)
appName <- 'ChroMo'

ui <-
  navbarPage(
    windowTitle = appName,
    title=div(style = "position: relative; top: -5px;", img(height = "35px", src='logo.png'), HTML(paste0('<b>',appName,'</b>'))),
    fluid = TRUE,
    collapsible = TRUE,
    theme = shinytheme("lumen"),
    tabPanel("Analysis", icon = icon("calculator"),
             sidebarLayout(
               sidebarPanel(
                 width = 3,
                 
                 # Upload
                 conditionalPanel(
                   condition = "input.tabs=='Data upload'",
                   h4("Data upload"),
                   h5("Upload from"),
                   radioButtons(
                     "data_input",
                     "Select a source:",
                     choices =
                       list(
                         "Local file" = 1,
                         "Web address" = 2,
                         "Biological Example" = 3,
                         "Synthetic Example (big!)" = 4
                       )
                     ,
                     selected =  1
                   ),
                   
                   conditionalPanel(
                     condition = "input.data_input=='1'",
                     fileInput(
                       "upload",
                       "Select a local file:",
                       accept = c(
                         'application/vnd.ms-excel',
                         'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                         '.xls',
                         '.xlsx',
                         '.csv',
                         '.tsv'
                       ),
                       multiple = TRUE
                     )
                   ),
                   conditionalPanel(condition = "input.data_input=='2'",
                                    textInput("URL", "URL", value = ""),
                                    NULL),
                   conditionalPanel(
                     condition = "output.fileUploaded == true",
                     h5("Time series meaning"),
                     selectizeInput(
                       "time_vars",
                       "Time variable",
                       choices = list("None" = "none"),
                       multiple = FALSE,
                       options = NULL
                     ),
                     selectizeInput(
                       "coord_vars",
                       "Coordinates variables",
                       choices = list("None" = "none"),
                       multiple = TRUE,
                       options = NULL
                     ),
                     numericInput(
                       "sampl_freq",
                       "Sampling period (s/frame)",
                       value = 60
                     ),
                     h5("Data blocking"),
                     selectizeInput(
                       "particle_vars",
                       "Particle identifier",
                       choices = list("None" = "none"),
                       multiple = FALSE,
                       options = NULL
                     ),
                     selectizeInput(
                       "grouping_vars",
                       "Group identifier",
                       choices = list("None" = "none"),
                       multiple = FALSE,
                       options = NULL
                     ),
                     checkboxInput(
                       inputId = "show_selected",
                       label = "Only show selected",
                       value = FALSE
                     ),
                     h5("Preprocess"),
                     radioButtons(
                       "norm_type",
                       "Data scaling:",
                       choices = list(
                         "None" = "none",
                         "Maximum-scaling" = "max",
                         "Minimum-scaling" = "min",
                         "0-1 scaling" = "zero_one",
                         "Z-score" = "z-score"
                       ),
                       selected = "none"
                     ),
                     conditionalPanel(
                      condition = "input.norm_type != 'none'",
                      selectizeInput(
                        "cols_normalize",
                        "Find columns to normalize",
                        choices = list("None" = "none"),
                        selected = "none",
                        multiple = TRUE,
                        options = NULL
                      )
                     ),
                     hr(),
                     checkboxInput(
                       inputId = "rotate_coords_seg",
                       label = "Rotate to maximum coordinates",
                       value = FALSE
                     ),
                     hr(),
                     checkboxInput(
                       inputId = "calc_veloc",
                       label = "Add vectorial and angular velocity",
                       value = FALSE
                     ),
                     checkboxInput(
                       inputId = "calc_more",
                       label = "Add cumulative distance",
                       value = FALSE
                     ),
                     hr(),
                     radioButtons(
                       "norm_time",
                       "Time normalization",
                       choices = list(
                         "None" = "none",
                         "Last is common" = "last",
                         "Last is common and zero" = "last_norm"
                       ),
                       selected = "none"
                     ),
                     hr(),
                     checkboxInput(
                       inputId = "exclude_cells",
                       label = "Exclude cells",
                       value = FALSE
                     ),
                     conditionalPanel(
                      condition = "input.exclude_cells == true && input.particle_vars != 'none'",
                      selectizeInput(
                        "cells_to_remove",
                        "Particles to remove",
                        choices = list("None" = "none"),
                        selected = "none",
                        multiple = TRUE,
                        options = NULL
                      ),
                     ),
                     hr(),
                     h4("Plot export"),
                          sliderInput(
                            "plot_width",
                            label = "Width (px)",
                            min = 100,
                            max = 1800,
                            value = c(600)
                          ),
                          sliderInput(
                              "plot_height",
                              label = "Height (px)",
                              min = 100,
                              max = 1800,
                              value = c(400)
                          )
                   ),
                   NULL
                 ),
                 
                 # Segments
                 conditionalPanel(
                   condition = "input.tabs=='Segments'",
                   
                   # Clusters: segmentation methodology
                   h4("Segments"),
                   h5("Segmentation criteria"),
                   selectizeInput(
                     "seg_method",
                     "Method",
                     choices = list("Lavielle" = "lavielle",
                                    "Multivariate" = "multiv"),
                     multiple = FALSE,
                     options = NULL
                   ),
                   selectizeInput(
                      "segment_variables",
                      "Columns to cluster",
                      choices = list("None" = "none"),
                      multiple = TRUE,
                      options = NULL
                   ),
                   actionButton('run_segmentation', 'Calculate segmentation'),
                   br(),
                   hr(),

                   # Clusters: what plots
                   conditionalPanel(
                     condition = "input.tabsegment == 'Discovery'",
                    sliderInput(
                      "seg_range",
                      label = "Cluster range",
                      min = 1,
                      max = 20,
                      value = c(1, 10)
                    ),
                    h6("Optimum cluster size will be selected"),


                    hr(),
                     h5("Test parameters"),
                     selectizeInput(
                       "signif_segment_test",
                       "Significance test",
                       choices = list(
                         "Likelihood Ratio" = "LRT",
                         "Chi-squared" = "Chisq",
                         "F-statistic" = "F",
                         "Rao-test" = "Rao"
                       ),
                       selected = "Likelihood",
                       multiple = FALSE,
                       options = NULL
                     ),
                     selectizeInput(
                       "signif_segment_adj",
                       "P-value adjustment",
                       choices = list(
                         "None" = "none",
                         "Bonferroni" = "fdr",
                         "Benjamini-Hochberg" = "holm"
                       ),
                       selected = "none",
                       multiple = FALSE,
                       options = NULL
                     )
                   ),
                   conditionalPanel(
                     condition = "input.tabsegment == 'Spectrum'",
                     h5("Spectrogram"),
                     sliderInput(
                       "spec_range",
                       label = "Spectrum comparison",
                       min = 0,
                       max = 1,
                       value = c(0.01, 0.6)
                     ),
                     checkboxInput(
                       inputId = "freq_spectrogram",
                       label = "Frequencies?",
                       value = FALSE
                     ),
                     checkboxInput(
                       inputId = "spectrum_individual",
                       label = "Individual?",
                       value = FALSE
                     ),
                     actionButton('run_spectrum', 'Calculate spectrum'),
                   ),
                   conditionalPanel(
                     condition = "input.tabsegment == 'Velocities'",
                     h5("Parameters"),
                     sliderInput(
                       "vel_ma",
                       label = "Smoothing",
                       min = 1,
                       max = 20,
                       value = 10
                     ),
                     checkboxInput(
                       inputId = "velocities_individual",
                       label = "Individual?",
                       value = FALSE
                     ),
                     actionButton('run_velocities', 'Calculate Velocities'),
                   ),
                   conditionalPanel(
                     condition = "input.tabsegment == 'Displacement'",
                     h5("Parameters"),
                     sliderInput(
                       "displ_k",
                       label = "Maximum lag (k)",
                       min = 1,
                       max = 40,
                       value = 20
                     ),
                     actionButton('run_msd', 'Calculate MSD'),
                   ),
                   conditionalPanel(
                     condition = "input.tabsegment == 'Individual'",
                     h5("What to show?"),
                     selectizeInput(
                       "cols_individual_plot",
                       "Values to plot",
                       choices = list("None" = "none"),
                       selected = "none",
                       multiple = TRUE,
                       options = NULL
                     ),
                     selectizeInput(
                       "part_individual_plot",
                       "Particles to plot",
                       choices = list("None" = "none"),
                       selected = "none",
                       multiple = TRUE,
                       options = NULL
                     ),
                     conditionalPanel(
                      condition = "input.part_individual_plot.length == 1",
                      checkboxInput(
                        inputId = "ind_specgram",
                        label = "Spectrogram",
                        value = FALSE
                      ),
                      conditionalPanel(
                        condition = "input.ind_specgram == true",
                        sliderInput(
                          "ind_spec_level",
                          label = "Significance",
                          min = 0,
                          max = 1,
                          value = 0.99
                        )
                      )
                     )
                   ),
                   NULL
                 ),
                 
                 # Motifs
                 conditionalPanel(
                   condition = "input.tabs=='Motifs'",
                   h4("Motifs"),
                   h5("Discovery parameters"),
                   selectizeInput(
                     "df_vars_motifs",
                     "Variable to analyze",
                     choices = list("None" = "none"),
                     selected = "none",
                     multiple = FALSE,
                     options = NULL
                   ),
                   conditionalPanel(
                     condition = "input.grouping_vars != 'none'",
                     selectizeInput(
                       "df_sequence_motifs",
                       "Analysis group",
                       choices = list("None" = "none"),
                       selected = "none",
                       multiple = TRUE,
                       options = NULL
                     ),
                     selectizeInput(
                       "df_query_motifs",
                       "Query group",
                       choices = list("None" = "none"),
                       selected = "none",
                       multiple = FALSE,
                       options = NULL
                     )
                   ),
                   sliderInput(
                     "motif_window",
                     label = "Window size",
                     min = 1,
                     max = 100,
                     value = 20
                   ),
                   sliderInput(
                     "motif_sample_pct",
                     label = "Samples to compute",
                     min = 0,
                     max = 1,
                     value = 1
                   ),
                   sliderInput(
                     "motif_correlation",
                     label = "Correlation",
                     min = 0,
                     max = 1,
                     value = 0.98
                   ),
                   sliderInput(
                     "motif_amount",
                     label = "Motifs to discover",
                     min = 1,
                     max = 6,
                     value = 3
                   ),
                   conditionalPanel(
                    condition = "input.tabmotifs=='Per segment'",
                    actionButton('run_motifs', 'Run Motif Discovery')
                   ),
                   hr(),
                   h5("Plots to show"),
                   checkboxInput(
                     inputId = "motif_show_mp",
                     label = "Matrix Profile",
                     value = FALSE
                   ),
                   checkboxInput(
                     inputId = "motif_show_motifs",
                     label = "Motifs",
                     value = TRUE
                   ),
                   checkboxInput(
                     inputId = "motif_show_discord",
                     label = "Discords",
                     value = FALSE
                   ),
                   NULL
                 ),
                 
                 # Correlations
                 conditionalPanel(
                   condition = "input.tabs=='Causality'",
                   h4("Causality analysis"),
                   selectizeInput(
                       "df_query_causality",
                       "Selection of variables",
                       choices = list("None" = "none"),
                       selected = "none",
                       multiple = TRUE,
                       options = NULL
                   ),
                   selectizeInput(
                       "df_sequence_causality",
                       "Data from group",
                       choices = list("None" = "none"),
                       selected = "none",
                       multiple = FALSE,
                       options = NULL
                     ),
                  # TODO: implement the causality per segment
                  #  checkboxInput(
                  #    inputId = "segment_correlate",
                  #    label = "Causality per segment",
                  #    value = FALSE
                  #  ),
                   hr(),
                   h4("Test parameters"),
                   sliderInput(
                      "pval_corr_range",
                      label = "Significance threshold (p-value)",
                      min = 0,
                      max = 0.25,
                      value = c(0.05)
                    ),
                   sliderInput(
                      "lags_corr_range",
                      label = "Maximum lags to explore",
                      min = 0,
                      max = 50,
                      value = c(10)
                    ),
                    sliderInput(
                      "perc_corr_range",
                      label = "Connection presence",
                      min = 0,
                      max = 1,
                      value = c(0.8)
                    ),
                    hr(),
                    actionButton('run_causal', 'Run causality analysis')
                 )
               ),
               
               mainPanel(
                 tags$head(tags$script(src = 'https://cdn.plot.ly/plotly-latest.min.js')),
                 tabsetPanel(
                   id = "tabs",
                   tabPanel("Data upload", icon = icon("upload"),
                      br(),
                      conditionalPanel(
                        condition = "output.fileUploaded == true",
                        fluidPage(
                          fluidRow(
                            column(4, h4("Uploaded data")),
                            column(1, downloadButton("downCSV", "Download data")),
                          )
                        ),
                        hr()
                      ),
                      dataTableOutput("data_uploaded")
                    ),
                   tabPanel(
                     "Segments",icon = icon("puzzle-piece"),
                     br(),
                     tabsetPanel(id = "tabsegment",
                      tabPanel("Discovery", icon = icon("microscope"),
                        h4("Group plot"),
                        downloadButton("downSegm", "Download figure"),
                        withSpinner(plotlyOutput("plot_segments", height = "auto")),
                        hr(),
                        h4("Clustered Segments significance"),
                        withSpinner(verbatimTextOutput('table_segments_significant'))
                      ),
                      tabPanel("Segment summary", icon = icon("book"),
                        h4("Summary table"),
                        withSpinner(dataTableOutput('table_segments'))
                      ),
                      tabPanel("Spectrum", icon = icon("wave-square"),
                        br(),
                        tabsetPanel(id = "spectrumtab",
                          tabPanel("Grouped plots", icon = icon("layer-group"),
                            h4("Spectral densities plot"),
                            downloadButton("downSpec", "Download figure"),
                            withSpinner(plotlyOutput("plot_spectrum", height = "auto")),
                            hr(),
                            h4("Spectral densities per group"),
                            downloadButton("downSpecSeg", "Download figure"),
                            withSpinner(plotlyOutput("plot_spectrum_segments", height = "auto"))
                          ),
                          tabPanel("Heatmaps", icon = icon("map"),
                            h4("Spectral heatmaps"),
                            downloadButton("downSpecH", "Download figure"),
                            withSpinner(plotOutput("plot_heatmap", height = "auto", width = 600)),
                            hr(),
                            h4("Spectral heatmaps per segment"),
                            downloadButton("downSpecSegH", "Download figure"),
                            withSpinner(plotOutput("plot_heatmap_segment", height = "auto", width = 600))
                          ),
                          tabPanel("Summary", icon = icon("book"),
                            h4("Spectral difference significance"),
                            withSpinner(verbatimTextOutput("summ_spectrum"))
                          )
                        )
                      ),
                      tabPanel("Velocities", icon = icon("tachometer-alt"),
                        br(),
                        tabsetPanel(id = "velocitytab",
                          tabPanel("Grouped plots",  icon = icon("layer-group"),
                            h4("Velocity densities"),
                            downloadButton("downVel", "Download figure"),
                            withSpinner(plotlyOutput("plot_velocities", height = "auto")),
                            hr(),
                            h4("Velocity per segment"),
                            downloadButton("downVelSeg", "Download figure"),
                            withSpinner(plotlyOutput("plot_velocities_segment", height = "auto"))
                          ),
                          tabPanel("Heatmaps",  icon = icon("map"),
                            h4("Velocity heatmaps"),
                            downloadButton("downVelH", "Download figure"),
                            withSpinner(plotOutput("plot_velocities_heat", height = "auto", width = 600)),
                            hr(),
                            h4("Velocity heatmaps per segment"),
                            downloadButton("downVelSegH", "Download figure"),
                            withSpinner(plotOutput("plot_velocities_heat_segment", height = "auto", width = 600))
                          ),
                          tabPanel("Summary",  icon = icon("book"),
                            h4("Velocities difference significance"),
                            withSpinner(verbatimTextOutput("summ_velocities"))
                          )
                        )
                      ),
                      tabPanel("Displacement", icon = icon("walking"),
                        h4("Mean Squared Displacement (MSD)"),
                        downloadButton("downMSD", "Download figure"),
                        withSpinner(plotlyOutput("plot_msd", height = "auto")),
                        hr(),
                        h4("MSD by segment"),
                        downloadButton("downMSDSeg", "Download figure"),
                        withSpinner(plotlyOutput("plot_msd_segment", height = "auto"))
                      ),
                      tabPanel("Individual", icon = icon("dice-one"),
                      h4("Individual plot"),
                      downloadButton("downIndiv", "Download figure"),
                      withSpinner(plotlyOutput("plot_individual", height = "auto")),
                      conditionalPanel(condition = "input.ind_specgram == true",
                                      fluidPage(
                                        fluidRow(
                                          column(4, h4("Individual Spectrogram")),
                                          column(1, downloadButton("downSpecgram", "Download figure")),
                                        )
                                      ),
                                      withSpinner(
                                        plotOutput("plot_morlet", height = "auto")
                                      ))
                      )
                     ),
                     NULL
                   ),
                   tabPanel(
                     "Motifs",icon = icon("hat-wizard"),
                     conditionalPanel(
                      condition = "input.df_vars_motifs != 'none'",
                      tabsetPanel(id = "tabmotifs",
                        tabPanel("Global", icon = icon("layer-group"),
                          br(),
                          conditionalPanel(
                            condition = "input.motif_show_mp == true",
                            fluidPage(
                              fluidRow(
                                column(4, h4("Global matrix profile")),
                                column(1, downloadButton("downMP", "Download figure")),
                              )
                            ),
                            withSpinner(plotOutput("plot_motifs_a", height = "auto")),
                            hr()
                          ),
                          conditionalPanel(
                            condition = "input.motif_show_motifs == true",
                            fluidPage(
                              fluidRow(
                                column(4, h4("Discovered motifs")),
                                column(1, downloadButton("downMotif", "Download figure")),
                              )
                            ),
                            withSpinner(plotOutput("plot_motifs_b", height = "auto")),
                            hr()
                          ),
                          conditionalPanel(
                            condition = "input.motif_show_discord == true",
                            fluidPage(
                              fluidRow(
                                column(4, h4("Discovered discords")),
                                column(1, downloadButton("downDiscord", "Download figure")),
                              )
                            ),
                            withSpinner(plotOutput("plot_motifs_c", height = "auto")),
                            hr()
                          )
                        ),
                        tabPanel("Per segment", icon = icon("dice-one"),
                          h4("Discovered per-segment motifs"),
                          h5("Only suitable-sized segments will be processed"),
                            withSpinner(uiOutput("tabs_motifs"))
                        )
                      )
                    )
                   ),
                   tabPanel(
                     "Causality", icon = icon("connectdevelop"),
                     br(),
                     tabsetPanel(id = "tabcorrs",
                        tabPanel("Global", icon = icon("layer-group"),
                          br(),
                          fluidPage(
                            fluidRow(
                              column(4, h4("PC algorithm")),
                              column(1, downloadButton("downPC", "Download figure")),
                            )
                          ),
                          withSpinner(forceNetworkOutput("plot_pc_corr", width = "600px", height = "600px")),
                          hr(),
                          fluidPage(
                            fluidRow(
                              column(4, h4("Variable Lag Transfer Entropy")),
                              column(1, downloadButton("downVLTE", "Download figure")),
                            )
                          ),
                          withSpinner(forceNetworkOutput("plot_vlte_corr", width = "600px", height = "600px")),
                          hr(),
                          fluidPage(
                            fluidRow(
                              column(4, h4("Variable Lag Granger Causality")),
                              column(1, downloadButton("downVLGC", "Download figure")),
                            )
                          ),
                          withSpinner(forceNetworkOutput("plot_vlgc_corr", width = "600px", height = "600px")),
                          hr(),
                          fluidPage(
                            fluidRow(
                              column(4, h4("Correlation"))
                              #column(1, downloadButton("downVLTE", "Download figure")),
                            )
                          ),
                          withSpinner(plotOutput("plot_cor", height = "auto"))
                        ),
                        tabPanel("Matrix", icon = icon("border-all"),
                          br(),
                          fluidPage(
                              fluidRow(
                                column(10, h4("Average adjacency matrices"))
                              )
                            ),
                          withSpinner(verbatimTextOutput("adj_matrix")),
                          hr()
                        )
                      )
                    )
                   )
                 )
               )
             ),
    tabPanel("About",icon = icon("info-circle"),
             sidebarLayout(
               sidebarPanel(
                 h4("About"),
                 "There are currently",
                 verbatimTextOutput("count"),
                 "session(s) connected to this app.",
               ),
               mainPanel(includeHTML("about.html"))
             ))
  )

server <- function(input, output, session) {
  isolate(vals$count <- vals$count + 1)
  vals$Datum <- FALSE
  
  observeEvent(input$change_scale, {
    if (input$change_scale == TRUE)  {
      updateCheckboxInput(session, "change_scale2", value = TRUE)
    } else if (input$change_scale == FALSE)   {
      updateCheckboxInput(session, "change_scale2", value = FALSE)
    }
  })
  
  observeEvent(input$change_scale2, {
    if (input$change_scale2 == TRUE)  {
      updateCheckboxInput(session, "change_scale", value = TRUE)
    } else if (input$change_scale2 == FALSE)   {
      updateCheckboxInput(session, "change_scale", value = FALSE)
    }
  })
  
  fileUploaded <- FALSE
  # Input data upload
  df_upload <- reactive({
    data <- NULL
    if (input$data_input == 1) {
      if (is.null(input$upload)) {
        # From uploaded file
        data <- NULL
      } else {
        isolate({
          # TODO: implement various csv files, each is a group
          if (input$upload$type == "text/csv") {
            # df_input <- read.csv(input$upload$datapath)
            df_input_list <- lapply(input$upload$datapath, read.csv)
          } else if (grepl("xls", input$upload$type) ||
                     grepl("Excel", input$upload$type)) {
            # df_input <- read_excel(input$upload$datapath)
            df_input_list <- lapply(input$upload$datapath, read_excel)
          }

          names(df_input_list) <- gsub(input$upload$name, pattern="\\..*", replacement="")

          df_input <- bind_rows(df_input_list, .id = "group.file.chromo")

          data <- df_input
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
      if (input$norm_type == "max") {
        norm_temp <-
          df_selected() %>% group_by(!!sym(input$particle_vars)) %>%
          mutate_at(cols_normalize, normalize.max)
      } else if (input$norm_type == "min") {
        norm_temp <-
          df_selected() %>% group_by(!!sym(input$particle_vars)) %>%
          mutate_at(cols_normalize, normalize.min)
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
    if (input$calc_veloc && length(coords) == 2) {
      norm_temp <- norm_temp %>% group_by(!!sym(input$particle_vars), !!sym(input$grouping_vars)) %>%
                    mutate(ang.speed =
                            angular_speed(data.frame(!!sym(coords[1]), !!sym(coords[2])),
                                          coord.names = coords),
                           cal.speed =
                            calculate.velocity(!!sym(coords[1]), !!sym(coords[2])),)
    }
    if (input$calc_more && length(coords) == 2) {
      norm_temp <- norm_temp %>% group_by(!!sym(input$particle_vars), !!sym(input$grouping_vars)) %>%
                    mutate(cum.dist =
                            calculate.cumdist(!!sym(coords[1]), !!sym(coords[2])))
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

    input_names_cols_calc <- c(
      "cols_individual_plot",
      "multiv_variables",
      "df_vars_motifs",
      "segment_variables",
      "df_query_causality"
    )

    for (u in input_names_cols) {
      updateSelectizeInput(session, u,
                           choices = var_list)
    }

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
      updateSelectizeInput(session,
                           "df_sequence_motifs",
                           choices = groups[, 1],
                           selected = groups[1, 1])
      updateSelectizeInput(session,
                           "df_sequence_causality",
                           choices = groups[, 1],
                           selected = groups[1, 1])
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


  segmentation_calc <- reactive({
    segments <- clusterize.segments(req(segmentation_raw()[[1]]),
                                    input$grouping_vars,
                                    input$particle_vars,
                                    segments = input$seg_range)
    segmentation_clusters(create.df.clusters(segments, segmentation_raw()[[2]]))
    return(segments)
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
        labs(x = "Time") +
        theme(legend.position = "none")
      
      #p <- ggplotly(p, autosize = TRUE,
      #      width = 600,
      #      height = 200*length(input$cols_individual_plot)*length(input$part_individual_plot))
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
      need(spectrum_global(), "\n\n\nPlease, run 'Calculate Spectrum' to continue.\nNo data is available for global spectrum"
      )
    )
    p <- spectral.plot(
      spectrum_global(),
      input$sampl_freq,
      input$freq_spectrogram,
      input$spectrum_individual
    )
    return(p) # in Plotly format
  })

  plot_spectrum_segments <- reactive({
    validate(
      need(spectrum_segments(), "\n\n\n\nPlease, run 'Calculate Spectrum' to continue.\nNo data is available segmented spectrum"
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
      need(spectrogram_global(), "\n\n\nPlease, run 'Calculate Spectrum' to continue."
      )
    )
    p <- heatmap.plot(
      spectrogram_global(),
      input$freq_spectrogram
    )
    return(p) # in Plotly format
  })

  plot_heatmap_segment <- reactive({
    validate(
      need(spectrogram_segments(), "\n\n\n\nPlease, run 'Calculate Spectrum' to continue."
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
      need(velocities_group(), "\n\n\nPlease, run 'Calculate Velocities' to continue."
      )
    )
    p <- heatmap.veloc.plot(
      velocities_group()
    )
    return(p) # in Plotly format
  })

  plot_velocities_heat_segment <- reactive({
    validate(
      need(velocities_segment(), "\n\n\nPlease, run 'Calculate Velocities' to continue."
      )
    )
    p <- heatmap.veloc.segmented.plot(
      velocities_segment()
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
  
  motifs_discovery <- reactive({
    if (input$df_vars_motifs == "none") {
      return(NULL)
    }
    
    df <- df_normalized()
    seqs <-
      df[df[[input$grouping_vars]] == input$df_sequence_motifs, input$df_vars_motifs]

    time <-
      df[df[[input$grouping_vars]] == input$df_sequence_motifs, input$time_vars]
    groups <-
      df[df[[input$grouping_vars]] == input$df_sequence_motifs, input$grouping_vars]
    
    motifs <-
      motif.discovery(
        seqs,
        win = input$motif_window,
        pct = input$motif_sample_pct,
        thr = input$motif_correlation,
        time = time,
        group = groups,
        nmotifs = input$motif_amount
      )
    return(motifs)
  })
  
  motifs_discovery_segmented <- reactiveVal()

  observeEvent(input$run_motifs, {
    req(segmentation_clusters())
    if (input$df_vars_motifs == "none") {
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

    segments.all <- segmentation_clusters()
    df <- df_process()

    progress <- Progress$new(session)
    progress$set(message = "Computing Causality", value = 0)

    # Temporally retrieving reactive values for future
    df_query_causality <- input$df_query_causality
    pval_corr_range <- input$pval_corr_range
    lags_corr_range <- input$lags_corr_range
    df_sequence_causality <- input$df_sequence_causality

    cor_group(calc.cor(df,
                          df_query_causality,
                          df_sequence_causality))

    causality.list <- causality.global(df,
                          df_query_causality,
                          pval_corr_range,
                          lags_corr_range,
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
    spec_signif <- spectrum.significance(
      spectrum_global(),
      input$spec_range
    )
    return(spec_signif)
  })

  plot_morlet <- reactive({
    plt_ml <- df_normalized()
    plt_ml <-
      plt_ml %>% filter(!!sym(input$particle_vars) == input$part_individual_plot)
    plt <-
      morlet.plot(plt_ml, input$cols_individual_plot[1], sig = input$ind_spec_level)
    return(plt)
  })

  velocities_group <- reactiveVal()
  velocities_segment <- reactiveVal()
  segmentation_clusters <- reactiveVal()

  observeEvent(input$run_velocities, {
    df <- df_process()
    df.clu <- segmentation_clusters()

    # Temporally retrieving reactive values for future
    particle_vars <- input$particle_vars
    coord_vars <- input$coord_vars
    vel_ma <- input$vel_ma

    future(seed=TRUE, {
      df_vt.smooth <- calc.velocities(df,
                                    coord_vars,
                                    vel_ma)
      df_vt.smooth

    }) %...>% velocities_group()

    progress <- Progress$new(session)
    progress$set(message = "Computing Velocities", value = 0)

    future(seed=TRUE, {

      updateProgress <- function(detail = NULL) {
        progress$inc(amount = 1/2, detail = detail)
      }

      df_seg_vels <- velocity.per.segment(df.clu,
                      coord_vars,
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

    msd_segment(calc.msd.segment(segments.all,
                                  coord_vars,
                                  displ_k))
  })
  
  plot_velocities <- reactive({
    validate(
      need(velocities_group(), "\n\n\nPlease, run 'Calculate Velocities' to continue."
      )
    )
    p <- plot.velocities(
      velocities_group(),
      input$velocities_individual
    )
    return(p)
  })

  plot_velocities_segment <- reactive({
    validate(
      need(velocities_segment(), "\n\n\nPlease, run 'Calculate Velocities' to continue."
      )
    )
    p <- plot.velocities.segment(
      velocities_segment()
    )
    return(p)
  })

  plot_msd <- reactive({
    validate(
      need(msd_group(), "\n\n\nPlease, run 'Calculate MSD' to continue."
      )
    )
    p <- plot.msd(
      msd_group()
    )
    return(p)
  })

  plot_msd_segment <- reactive({
    validate(
      need(msd_segment(), "\n\n\n\nPlease, run 'Calculate MSD' to continue."
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
    significant.velocities(velocities_group())
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

  output$summ_velocities <- renderPrint({
    velocities_significance()
  })

  output$tabs_motifs <- renderUI({
    seg.toplot <- req(plot_motifs_segmented())
    range <- c(1, 180)
    for (i in 1:length(seg.toplot)) {
      local({
        my_i <- i
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
  
}

shinyApp(ui = ui, server = server)