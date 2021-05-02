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
                   )