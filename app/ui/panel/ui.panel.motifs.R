chromo.ui.panel.motifs <- tabPanel(
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
                        tabPanel("Per cluster", icon = icon("dice-one"),
                          h4("Discovered per-cluster motifs"),
                          h5("Only suitable-sized clusters will be processed"),
                            withSpinner(uiOutput("tabs_motifs"))
                        )
                      )
                    )
                   )