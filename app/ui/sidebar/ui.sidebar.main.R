source("ui/sidebar/ui.sidebar.upload.R")
source("ui/sidebar/ui.sidebar.segment.R")
source("ui/sidebar/ui.sidebar.motifs.R")
source("ui/sidebar/ui.sidebar.causality.R")

chromo.sidebar <- sidebarPanel(
                 width = 3,
                 
                 chromo.ui.sidebar.upload,
                 chromo.ui.sidebar.segment,
                 chromo.ui.sidebar.motifs,
                 chromo.ui.sidebar.causality
               )