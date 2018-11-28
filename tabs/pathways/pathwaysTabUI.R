tabs.pathways.moduleUI <- function (id) {

  ns <- NS(id)

  fluidRow(
    box(
      title = "Open in Reactome browser",
      solidHeader = TRUE,
      p("You can perform a pathway overrepresentation analysis here. Gene symbols are uploaded to Reactome's servers and the analysis website is opened in a new tab."),
      hr(),
      actionButton(ns("openReactomeButton"), label = "Open in  Reactome", class="btn-info btn-block"),
      uiOutput(ns("reactomeOpeningScript"))
    )
  )

}
