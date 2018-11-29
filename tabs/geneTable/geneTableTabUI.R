tabs.geneTable.moduleUI <- function (id) {

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
        util.tableDownload.moduleUI(ns("geneDownload"), "genes"),
        footer = "This download contains the filtered genes only. Please note that only the sidebar filters apply to this file, not the sorting orders or search bar."
      ),
      box(
        width = 3,
        solidHeader = TRUE,
        title = "Copy gene names to clipboard",
        uiOutput(ns("geneNamesToClipboard")),
        footer = "You can use this button to compare the genes names from this table to other variant tables using the gene comparison tab."
      )
    )
  )
}
