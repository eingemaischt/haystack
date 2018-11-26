tableDownloadUI <- function (id, title) {

  ns <- NS(id)

  return(tagList(
    tags$div(
      downloadButton(ns("csvDownload"), label = paste0("Download ", title, " as CSV"))
    ),
    tags$div(
      downloadButton(ns("xlsxDownload"), label = paste0("Download ", title, " as XLSX"))
    )
  ))

}