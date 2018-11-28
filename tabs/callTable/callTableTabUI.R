tabs.callTable.moduleUI <- function (id) {

  ns <- NS(id)

  tagList(
    fluidRow(
      box(
        style = "overflow-x: scroll",
        title = "Variant Call Table",
        width = 12,
        solidHeader = TRUE,
        collapsible = TRUE,
        DT::dataTableOutput(ns("callTable"))
      )
    ),
    fluidRow(
      box(
        title = "Upload",
        width = 6,
        p("You can upload your variant call table here. Make sure to upload only correctly formatted files (comma separated, periods as decimal seperators, variant table as first sheet for .xlsx files). Both .csv and .xlsx files are supported."),
        hr(),
        fileInput(ns("callFile"), "Choose file",
                  accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv",
                    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    ".xlsx")
        )
      ),
      box(
        title = "Download",
        width = 6,
        style = "height: 100%",
        p("You can download the filtered variant call table here. Please note that only the sidebar filters apply to this file, not the sorting orders or search bar."),
        hr(),
        util.tableDownload.moduleUI(ns("filteredCallDownload"), "filtered calls")
      )
    ),
    fluidRow(
      box(
        title = "Select displayed columns",
        p("Select which columns to display here. Beware to not disable columns that you want to use otherwise."),
        checkboxGroupInput(ns("selectedColumns"), label = "Show columns: ", inline = TRUE)
      ),
      box(
        width = 3,
        title = "Symbols with no expression",
        footer = "These symbols have a valid HGNC symbol but show no expression data in any tissue (or all tissues show tpm values of 0 / 'not detected' protein level).",
        tags$div(
          style = "overflow-y: scroll; max-height: 400px",
          verbatimTextOutput(ns("unexpressedSymbols"))
        ),
        tags$div(
          "Total: ",
          textOutput(ns("totalUnexpressedSymbols"), inline = TRUE)
        ),
        uiOutput(ns("unexpressedToClipboard"))
      ),
      box(
        width = 3,
        title = "Unrecognized symbols",
        footer = "These symbols are no valid (up-to-date) HGNC symbols. Also, no valid previous or alias name could be found.",
        tags$div(
          style = "overflow-y: scroll; max-height: 400px",
          verbatimTextOutput(ns("unrecognizedSymbols"))
        ),
        tags$div(
          "Total: ",
          textOutput(ns("totalUnrecognizedSymbols"), inline = TRUE)
        ),
        uiOutput(ns("unrecognizedToClipboard"))
      )
    )
  )

}
