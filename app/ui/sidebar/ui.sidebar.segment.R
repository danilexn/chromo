chromo.ui.sidebar.segment <- function(request) {
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
                     sliderInput(
                       "spec_ticks",
                       label = "Number of ticks (x-axis)",
                       min = 0,
                       max = 20,
                       value = 5
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
                     condition = "input.tabsegment == 'Distributions'",
                     h4("Distributions"),
                     selectizeInput(
                      "density_variables",
                      "Input column",
                      choices = list("None" = "none"),
                      multiple = FALSE,
                      options = NULL
                   ),
                     h5("Parameters"),
                     sliderInput(
                       "vel_ma",
                       label = "Smoothing",
                       min = 1,
                       max = 20,
                       value = 10
                     ),
                     checkboxInput(
                       inputId = "densities_individual",
                       label = "Individual?",
                       value = FALSE
                     ),
                     actionButton('run_densities', 'Calculate Density'),
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
                 )
}