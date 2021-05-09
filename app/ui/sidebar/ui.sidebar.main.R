source("ui/sidebar/ui.sidebar.upload.R")
source("ui/sidebar/ui.sidebar.segment.R")
source("ui/sidebar/ui.sidebar.motifs.R")
source("ui/sidebar/ui.sidebar.causality.R")

chromo.sidebar <- function(request) {
              sidebarPanel(
                 width = 3,
                 h4("Save analysis"),
                 bookmarkButton(),
                 hr(),
                 chromo.ui.sidebar.upload(request),
                 chromo.ui.sidebar.segment(request),
                 chromo.ui.sidebar.motifs(request),
                 chromo.ui.sidebar.causality(request)
               )
}