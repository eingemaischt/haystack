util.tableDownload.moduleUI <- function (id, title, showXslx = TRUE) {

  ns <- NS(id)

  csvDiv <- tags$div(
    downloadButton(ns("csvDownload"), label = paste0("Download ", title, " as CSV"))
  )

  xlsxDiv <- tags$div(
    downloadButton(ns("xlsxDownload"), label = paste0("Download ", title, " as XLSX"))
  )

  if (!showXslx) return(csvDiv)

  return(tagList(csvDiv, xlsxDiv))
}
