source("ui/panel/ui.panel.upload.R")
source("ui/panel/ui.panel.segment.R")
source("ui/panel/ui.panel.motifs.R")
source("ui/panel/ui.panel.causality.R")

chromo.mainpanel <- mainPanel(
                 tags$head(tags$script(src = 'https://cdn.plot.ly/plotly-latest.min.js')),
                 tabsetPanel(
                   id = "tabs",
                   chromo.ui.panel.upload,
                   chromo.ui.panel.segment,
                   chromo.ui.panel.motifs,
                   chromo.ui.panel.causality
                   )
                 )