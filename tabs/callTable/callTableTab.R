tabs.callTable.isValidTable <- function (dt) {

  characterColumns <- c("Sample", "Studie", "Chr", "Symbol", "HGVSc", "Consequence", "Genotype")
  numericColumns <- c("Position", "AF Popmax", "Read depth", "Variant depth")
  expectedColumns <- c(characterColumns, numericColumns)

  if (any(!expectedColumns %in% colnames(dt))) return(FALSE)

  actualCharacterColumnTypes <- sapply(dt[,..characterColumns], class)
  actualNumericColumnTypes <- sapply(dt[,..numericColumns], class)

  return(
    all(actualCharacterColumnTypes == "character") &&
      all(actualNumericColumnTypes %in% c("numeric", "integer", "logical"))
  )

}

tabs.callTable.readCallTable <- function (fileName) {

  lowercaseFileName <- tolower(fileName)

  if (endsWith(lowercaseFileName, ".csv")) {

    ct <- fread(fileName)

  } else if (endsWith(lowercaseFileName, ".xlsx")) {

    # xlsx files are just parsed as character data frame,
    # so we need to convert column types to data tables.
    # the simplest way is to write the .xlsx to .csv and
    # then read using fread
    xlsxTable <- read.xlsx(fileName, check.names = FALSE)

    # check.names is still replacing spaces with dots. We replace them
    # to have consistent column names, see also:
    # https://github.com/awalker89/openxlsx/issues/102
    colnames(xlsxTable) <- gsub("\\.", " ", colnames(xlsxTable))

    tmpFile <- tempfile()
    fwrite(xlsxTable, tmpFile)
    ct <- fread(tmpFile)
    file.remove(tmpFile)

  } else {
    return(NULL)
  }

  if (is.numeric(ct$Chr)) {
    ct$Chr <- as.character(ct$Chr)
  }

  if (tabs.callTable.isValidTable(ct)) {
    return(ct)
  } else {
    return(NULL)
  }

}

tabs.callTable.module <- function (input, output, session) {

  fullCallTable <- reactiveVal()
  filteredCallTable <- reactiveVal()
  expressionFilter <- reactiveVal()

  callModule(tableDownload, "filteredCallDownload", filteredCallTable, "filtered-calls-")

  observeEvent(input$callFile, {

    req(input$callFile)

    progress <- shiny::Progress$new()
    progress$set(message = "Uploading table", value = .5)

    ct <- tabs.callTable.readCallTable(input$callFile$datapath)

    if(is.null(ct)) {
      util.showErrorModal("Could not correctly read table. Make sure there are no trailing rows, all relevant columns are present and dots are used as decimal seperators.")
      progress$close()
    }

    req(ct)

    fullCallTable(ct)

    updateCheckboxGroupInput(session, "selectedColumns", choices = colnames(ct), selected = colnames(ct))

    recognizedSymbolIndices <- annotation.symbolToIndexMap[[ct$Symbol]]
    symbolsAreRecognized <- !is.na(recognizedSymbolIndices)
    recognizedSymbols <- unique(annotation.geneTable$symbol[recognizedSymbolIndices[symbolsAreRecognized]])

    unexpressedSymbols <- recognizedSymbols[!recognizedSymbols %in% annotation.gtexExpression$symbol]
    unrecognizedSymbols <- unique(ct$Symbol[!symbolsAreRecognized])

    unrecognizedSymbolsText <- paste0(unrecognizedSymbols, collapse = "\n")
    unexpressedSymbolsText <- paste0(unexpressedSymbols, collapse = "\n")

    symbolsInTotal <- length(unique(ct$Symbol))
    unrecognizedInTotal <- length(unrecognizedSymbols)
    unexpressedInTotal <- length(unexpressedSymbols)

    fractionUnrecognized <- round(unrecognizedInTotal / symbolsInTotal, digits = 4) * 100
    fractionUnexpressed <- round(unexpressedInTotal / symbolsInTotal, digits = 4) * 100

    output$unrecognizedSymbols <- renderText(unrecognizedSymbolsText)
    output$unexpressedSymbols <- renderText(unexpressedSymbolsText)

    output$totalUnrecognizedSymbols <- renderText(paste0(unrecognizedInTotal, " (", fractionUnrecognized, "%)"))
    output$totalUnexpressedSymbols <- renderText(paste0(unexpressedInTotal, " (", fractionUnexpressed, "%)"))

    output$unrecognizedToClipboard <-util.copyToClipboardButton(unrecognizedSymbolsText, "unrecognizedbtn")
    output$unexpressedToClipboard <-util.copyToClipboardButton(unexpressedSymbolsText, "unexpressedbtn")

    progress$close()
  })


  output$callTable <- DT::renderDataTable({

    req(filteredCallTable())
    return(DT::datatable(filteredCallTable(),
                         selection = "single",
                         style = "bootstrap",
                         class = DT:::DT2BSClass(c("compact", "hover", "stripe")),
                         options = list(columnDefs = list(list(
                           targets = seq_len(ncol(filteredCallTable())),
                           render = DT::JS(
                             "function(data, type, row, meta) {",
                             "if (data == null) return '';",
                             "return type === 'display' && data.length > 20 ?",
                             "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                             "}"
                           )))
                         )))
  })

  observeEvent(input$callTable_rows_selected,
                util.detailModal.showDetailModal(
                 filteredCallTable()[input$callTable_rows_selected]$Symbol,
                 filteredCallTable,
                 expressionFilter(),
                 output,
                 "callTableTab")
  )

  observeEvent(input$modalOkBtn, {
    removeModal()
  })

  return(list(
    filteredCallTable = filteredCallTable,
    fullCallTable = fullCallTable,
    selectedColumns = reactive({ input$selectedColumns }),
    expressionFilter = expressionFilter
  ))
}
