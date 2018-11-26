geneTableTabUI <- function (id) {

  ns <- NS(id)

  tagList(
    fluidRow(
      box(
        style = "overflow-x: scroll",
        title = "Gene Table",
        width = 12,
        solidHeader = TRUE,
        collapsible = TRUE,
        DT::dataTableOutput(ns("geneTable"))
      )
    ),
    fluidRow(
      box(
        width = 6,
        solidHeader = TRUE,
        title = "Download",
        tableDownloadUI(ns("geneDownload"), "genes"),
        footer = "This download contains the filtered genes only. Please note that only the sidebar filters apply to this file, not the sorting orders or search bar."
      )
    )
  )
}