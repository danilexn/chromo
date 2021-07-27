chromo.ui.sidebar.motifs <- function(request) {
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
                    condition = "input.tabmotifs=='Per cluster'",
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
                 )
}