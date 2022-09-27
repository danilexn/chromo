# From https://stackoverflow.com/questions/46776538/r-shiny-combine-the-slider-bar-with-a-text-input-to-make-the-slider-bar-more-u

controlledSliderUI <- function(id, label){
  ns = NS(id)
  wellPanel(
    sliderInput(ns("slider"), label = label, 0, 1, c(0, 1)),
    textInput(ns("min"), "min", 0, "50%"),
    textInput(ns("max"), "max", 100, "50%"),
    actionButton(ns("update"), "update slider")
  )
}

controlledSlider <- function(input, output, session, min, max, value){
  reactiveRange <- reactiveValues(min = value[1], max = value[2])
  updateSliderInput(session, "slider", min = min, max = max)
  
  ## observe slider
  observeEvent(input$slider,{
    reactiveRange$min <- input$slider[1]
    reactiveRange$max <- input$slider[2]
  }, ignoreInit = TRUE)
  
  ## observe button
  observeEvent(input$update,{reactiveRange$min <- as.numeric(input$min)})
  observeEvent(input$update,{reactiveRange$max <- as.numeric(input$max)})
  
  ## observe reactive
  observeEvent({reactiveRange$min; reactiveRange$max},{
    updateSliderInput(
      session, "slider", value = c(reactiveRange$min, reactiveRange$max))
    updateTextInput(session, "min", value = reactiveRange$min)
    updateTextInput(session, "max", value = reactiveRange$max)
  })
  
  return(reactiveRange)
}