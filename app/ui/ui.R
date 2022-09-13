source("ui/sidebar/ui.sidebar.main.R")
source("ui/panel/ui.panel.main.R")
source("ui/custom_ui.R")

appName <- 'ChroMo'

ui <- function(request) {
  navbarPage(
    windowTitle = appName,
    title=div(style = "position: relative; top: -5px;", img(height = "35px", src='logo.png'), HTML(paste0('<b>',appName,'</b>'))),
    fluid = TRUE,
    collapsible = TRUE,
    theme = shinytheme("lumen"),
    tabPanel("Analysis", icon = icon("calculator"),
             sidebarLayout(
               chromo.sidebar(request),
               chromo.mainpanel
               )
             ),
    tabPanel("About",icon = icon("info-circle"),
             sidebarLayout(
               sidebarPanel(
                 h4("About"),
                 "There are currently",
                 verbatimTextOutput("count"),
                 "session(s) connected to this app.",
               ),
               mainPanel(includeHTML("about.html"))
             ))
  )
}