csvWriteHandler <- function (dataTable, fileName) {

  fwrite(dataTable, fileName)

}

xlsxWriteHandler <- function (dataTable, fileName) {

  write.xlsx(dataTable, fileName, asTable = TRUE)

}

handleTableDownload <- function (tableReactiveValue, filePrefix, fileExtension, writeHandler) {

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

tableDownload <- function (input, output, session, tableReactiveValue, filePrefix) {

  output$csvDownload <- handleTableDownload(tableReactiveValue, filePrefix, ".csv", csvWriteHandler)
  output$xlsxDownload <- handleTableDownload(tableReactiveValue, filePrefix, ".xlsx", xlsxWriteHandler)

}