annotationsTabUI <- function (id) {

  ns <- NS(id)

  tagList(
    fluidRow(
      box(
        style = "overflow-x: scroll",
        title = "Annotation Table",
        width = 12,
        solidHeader = TRUE,
        footer = "This window shows the annotation data used for this app, it is extracted from the HGNC.",
        DT::dataTableOutput(ns("annotationTable"))
      )
    ),
    fluidRow(
      box(
        width = 6,
        solidHeader = TRUE,
        footer = "Please note that downloads always contain the whole tables. Filters and sorting orders do not apply.",
        title = "Download",
        tableDownloadUI(ns("annotationDownload"), "gene annotations"),
        hr(),
        tableDownloadUI(ns("gtexDownload"), "GTEx table"),
        hr(),
        tableDownloadUI(ns("hpaRnaDownload"), "HPA RNA table"),
        hr(),
        tableDownloadUI(ns("hpaProteinDownload"), "HPA protein table")
      )
    )
  )

}