chromo.ui.panel.causality <- tabPanel(
                     "Causality", icon = icon("connectdevelop"),
                     br(),
                     tabsetPanel(id = "tabcorrs",
                        tabPanel("PC-alg", icon = icon("project-diagram"),
                          br(),
                          fluidPage(
                            fluidRow(
                              column(4, h4("PC algorithm")),
                              column(1, downloadButton("downPC", "Download figure")),
                            )
                          ),
                          withSpinner(forceNetworkOutput("plot_pc_corr", width = "600px", height = "600px")),
                        ),
                        tabPanel("VLTE", icon = icon("project-diagram"),
                          br(),
                          fluidPage(
                            fluidRow(
                              column(4, h4("Variable Lag Transfer Entropy")),
                              column(1, downloadButton("downVLTE", "Download figure")),
                            )
                          ),
                          withSpinner(forceNetworkOutput("plot_vlte_corr", width = "600px", height = "600px")),
                        ),
                        tabPanel("VLGC", icon = icon("project-diagram"),
                          br(),
                          fluidPage(
                            fluidRow(
                              column(4, h4("Variable Lag Granger Causality")),
                              column(1, downloadButton("downVLGC", "Download figure")),
                            )
                          ),
                          withSpinner(forceNetworkOutput("plot_vlgc_corr", width = "600px", height = "600px")),
                        ),
                        tabPanel("Correlation", icon = icon("chart-area"),
                          br(),
                          fluidPage(
                            fluidRow(
                              column(4, h4("Correlation"))
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