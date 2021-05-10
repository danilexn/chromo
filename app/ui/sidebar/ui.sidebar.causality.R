chromo.ui.sidebar.causality <- function(request) {
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
}