chromo.ui.panel.upload <- tabPanel("Data upload", icon = icon("upload"),
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
                    )