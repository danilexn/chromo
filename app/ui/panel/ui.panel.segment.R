chromo.ui.panel.segment <- tabPanel(
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
                      tabPanel("Cluster summary", icon = icon("book"),
                        h4("Summary table"),
                        h5("Displaying z-scaled means"),
                        withSpinner(dataTableOutput('table_segments')),
                        h4("Summary table per group"),
                        h5("Displaying z-scaled means"),
                        withSpinner(dataTableOutput('table_segments_grouped'))
                      ),
                      tabPanel("Spectrum", icon = icon("wave-square"),
                        br(),
                        tabsetPanel(id = "spectrumtab",
                          tabPanel("Grouped plots", icon = icon("layer-group"),
                            h4("Spectral densities plot"),
                            fluidRow(
                              column(4, downloadButton("downSpec", "Download figure"),),
                              column(4, downloadButton("downSpecTable", "Download CSV"),),
                            ),
                            withSpinner(plotlyOutput("plot_spectrum", height = "auto")),
                            hr(),
                            h4("Spectral densities per group"),
                            fluidRow(
                              column(4, downloadButton("downSpecSeg", "Download figure"),),
                              column(4, downloadButton("downSpecSegTable", "Download CSV"),),
                            ),
                            withSpinner(plotlyOutput("plot_spectrum_segments", height = "auto"))
                          ),
                          tabPanel("Heatmaps", icon = icon("map"),
                            h4("Spectral heatmaps"),
                            downloadButton("downSpecH", "Download figure"),
                            withSpinner(plotOutput("plot_heatmap", height = "auto", width = 600)),
                            hr(),
                            h4("Spectral heatmaps per cluster"),
                            downloadButton("downSpecSegH", "Download figure"),
                            withSpinner(plotOutput("plot_heatmap_segment", height = "auto", width = 600))
                          ),
                          tabPanel("Summary", icon = icon("book"),
                            h4("Spectral difference significance"),
                            h5("Global"),
                            withSpinner(verbatimTextOutput("summ_spectrum")),
                            h5("Per segment"),
                            withSpinner(verbatimTextOutput("summ_spectrum_segment"))
                          )
                        )
                      ),
                      tabPanel("Distributions", icon = icon("tachometer-alt"),
                        br(),
                        tabsetPanel(id = "densitytab",
                          tabPanel("Grouped plots",  icon = icon("layer-group"),
                            h4("Global distribution"),
                            fluidRow(
                              column(4, downloadButton("downVel", "Download figure"),),
                              column(4, downloadButton("downVelTable", "Download CSV"),),
                            ),
                            withSpinner(plotlyOutput("plot_velocities", height = "auto")),
                            hr(),
                            h4("Density per cluster"),
                            fluidRow(
                              column(4, downloadButton("downVelSeg", "Download figure"),),
                              column(4, downloadButton("downVelSegTable", "Download CSV"),),
                            ),
                            withSpinner(plotlyOutput("plot_velocities_segment", height = "auto"))
                          ),
                          tabPanel("Heatmaps",  icon = icon("map"),
                            h4("Global heatmaps"),
                            downloadButton("downVelH", "Download figure"),
                            withSpinner(plotOutput("plot_velocities_heat", height = "auto", width = 600)),
                            hr(),
                            h4("Heatmaps per cluster"),
                            downloadButton("downVelSegH", "Download figure"),
                            withSpinner(plotOutput("plot_velocities_heat_segment", height = "auto", width = 600))
                          ),
                          tabPanel("Summary",  icon = icon("book"),
                            h4("Density difference significance"),
                            h5("Global"),
                            withSpinner(verbatimTextOutput("summ_velocities")),
                            h5("Per segment"),
                            withSpinner(verbatimTextOutput("summ_velocities_segment"))
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
                   )