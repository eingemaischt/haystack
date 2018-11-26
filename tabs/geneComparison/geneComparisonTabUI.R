geneComparisonTabUI <- function (id) {

  ns <- NS(id)

  tagList(
    fluidRow(
      box(
        width = 12,
        title = "Venn Diagram",
        footer = "You can save this diagram using right click -> 'Save image as'.",
        plotOutput(ns("geneComparisonVenn"))
      )
    ),
    fluidRow(
      box(
        width = 4,
        title = "Symbols only in our data",
        tags$div(
          style = "overflow-y: scroll; max-height: 400px",
          verbatimTextOutput(ns("ourSymbols"))
        ),
        hr(),
        uiOutput(ns("ourSymbolsToClipboard"))
      ),
      box(
        width = 4,
        title = "Intersecting symbols",
        tags$div(
          style = "overflow-y: scroll; max-height: 400px",
          verbatimTextOutput(ns("intersectingSymbols"))
        ),
        hr(),
        uiOutput(ns("intersectingSymbolsToClipboard"))
      ),
      box(
        width = 4,
        title = "Symbols only in their data",
        tags$div(
          style = "overflow-y: scroll; max-height: 400px",
          verbatimTextOutput(ns("theirSymbols"))
        ),
        hr(),
        uiOutput(ns("theirSymbolsToClipboard"))
      )
    ),
    fluidRow(
      box(
        width = 4,
        title = "Gene List Upload",
        footer = "You can upload a gene list here to compare to the genes found in the variant table. Every line must contain exactly one HGNC symbol.",
        fileInput(ns("geneComparisonListUpload"), label = "Upload genes for comparison")
      ),
      box(
        width = 3,
        title = "Filters",
        footer = "Choose here if sidebar filters should be applied to symbols in our data or if all symbols should be considered regardless of filtering.",
        checkboxInput(ns("geneComparisonUseFilters"), value = TRUE, label = "Apply sidebar filters to our symbols")
      )
    )
  )

}