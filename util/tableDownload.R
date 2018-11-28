util.tableDownload.csvWriteHandler <- function (dataTable, fileName) {

  fwrite(dataTable, fileName)

}

util.tableDownload.xlsxWriteHandler <- function (dataTable, fileName) {

  write.xlsx(dataTable, fileName, asTable = TRUE)

}

util.tableDownload.handleTableDownload <- function (tableReactiveValue, filePrefix, fileExtension, writeHandler) {

  return(downloadHandler(
    filename = function() {
      paste0(filePrefix, Sys.time(), fileExtension)
    },
    content = function (con) {

      if(is.null(tableReactiveValue())) {

        writtenData <- data.table("NO DATA AVAILABLE")

      } else {

        writtenData <- tableReactiveValue()

      }

      writeHandler(writtenData, con)
    }
  ))

}

util.tableDownload.module <- function (input, output, session, tableReactiveValue, filePrefix) {

  output$csvDownload <- util.tableDownload.handleTableDownload(tableReactiveValue, filePrefix, ".csv", util.tableDownload.csvWriteHandler)
  output$xlsxDownload <- util.tableDownload.handleTableDownload(tableReactiveValue, filePrefix, ".xlsx", util.tableDownload.xlsxWriteHandler)

}
