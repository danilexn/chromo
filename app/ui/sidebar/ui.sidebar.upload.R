chromo.ui.sidebar.upload <- function(request) {
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
                         "Last is common and zero" = "last",
                         "Last is common" = "last_norm"
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
                   )
                 )
}