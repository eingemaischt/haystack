options(shiny.maxRequestSize=100*1024^2)

shiny.huge.detailDivElement <- function (label, value) {

  return(tags$div(
    tags$b(label),
    tags$em(value)
  ))

}

shiny.huge.modalExpressionPlot <- function (expressions, expressionFilter, title) {

  ### WORKAROUND
  ### when no expressions are found, just render some dummy plot so no errors are thrown
  ### see https://stackoverflow.com/questions/19918985/r-plot-only-text
  if (nrow(expressions) == 0) {

    errorMessage <- "No information available"

    return(renderPlot(
      ggplot() + annotate("text", x = 4, y = 25, size = 8, label = errorMessage)+
        theme_bw() +
        theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
    ))
  }
  ### END WORKAROUND

  return(renderPlot(
    ggplot(expressions, aes(x = reorder(tissue, value, FUN = max), y = value, fill = tissue %in% expressionFilter, colour = tissue %in% expressionFilter)) +
      geom_col() +
      coord_flip() +
      ggtitle(title) +
      xlab("Tissue") +
      ylab("Expression value") +
      scale_fill_manual(values = c("TRUE" = "springgreen", "FALSE" = "steelblue")) +
      scale_colour_manual(values = c("TRUE" = "springgreen4", "FALSE" = "steelblue4")) +
      guides(fill = FALSE, colour = FALSE),
    res = 90
    )
  )

}

shiny.huge.geneExpressionModal <- function (selectedSymbol, callTableReactiveVal, input, output) {

  return({

    req(callTableReactiveVal())

    matchingGene <- shiny.huge.geneTable[shiny.huge.symbolToIndexMap[[selectedSymbol]]]

    showModal(modalDialog(
      shiny.huge.detailDivElement("Location:", matchingGene$location),
      shiny.huge.detailDivElement("Gene (found in table):", selectedSymbol),
      shiny.huge.detailDivElement("Gene (current HGNC symbol):", matchingGene$symbol),
      shiny.huge.detailDivElement("Gene full name:", matchingGene$name),
      shiny.huge.detailDivElement("Gene family:", matchingGene$gene_family),
      shiny.huge.detailDivElement("Gene description:", matchingGene$description),
      shiny.huge.detailDivElement("Location type:", matchingGene$locus_type),
      shiny.huge.detailDivElement("Ensembl ID:", matchingGene$ensembl_gene_id),
      tags$hr(),
      tags$b("Expression:"),
      tabsetPanel(
        tabPanel("GTex (TPM scaled)", plotOutput("modalGTExScaledExpression", height = "640px")),
        tabPanel("HPA RNA (TPM scaled)", plotOutput("modalHpaRnaScaledExpression", height = "640px")),
        tabPanel("HPA Protein", plotOutput("modalHpaProteinExpression", height = "640px")),
        tabPanel("Mutation types", tags$div(
          style = "overflow-x: scroll",
          tags$h5(paste0("Mutations for ", selectedSymbol), ":"),
          tableOutput("modalMutationTypes")
          )
        )
      ),
      title = "Details",
      footer = actionButton("modalOkBtn", label = "OK", icon = icon("ok")),
      size = "l",
      easyClose = TRUE
    ))

    callsInGene <- callTableReactiveVal()[Symbol == selectedSymbol]

    uniqueMutations <- callsInGene[,
      list("Samples" = length(unique(Sample))),
      by = list(HGVSc, Chr, Position, Consequence, `AF Popmax`)
    ]

    matchingGtexExpression <- shiny.huge.gtexExpression[symbol == matchingGene$symbol]
    matchingHpaRnaExpression <- shiny.huge.hpaRnaExpression[symbol == matchingGene$symbol]
    matchingHpaProteinExpression <- shiny.huge.hpaProteinExpession[symbol == matchingGene$symbol]

    scaledGtexValues <- matchingGtexExpression[,list(tissue = tissue, value = tpm_scaled)]

    scaledHpaRnaValues <- matchingHpaRnaExpression[,list(tissue = tissue, value = tpm_scaled)]

    hpaProteinValues <- matchingHpaProteinExpression[,list(value = max(level)),by = tissue]

    output$modalGTExScaledExpression <- shiny.huge.modalExpressionPlot(scaledGtexValues, input$expressions, paste(selectedSymbol, "GTEx data (scaled TPM)", sep = ": "))
    output$modalHpaRnaScaledExpression <- shiny.huge.modalExpressionPlot(scaledHpaRnaValues, input$expressions, paste(selectedSymbol, "HPA RNA data (scaled TPM)", sep = ":"))
    output$modalHpaProteinExpression <- shiny.huge.modalExpressionPlot(hpaProteinValues, input$expressions, paste(selectedSymbol, "HPA Protein data (levels)", sep = ":"))
    output$modalMutationTypes <- renderTable(uniqueMutations, spacing = "xs")

  })

}

shiny.huge.showErrorModal <- function (errorMessage, session) {

  showModal(modalDialog(
    tags$div(errorMessage),
    title = "ERROR",
    easyClose = TRUE
  ))

}



shiny.huge.copyToClipboardButton <- function (text, id) {
  return(renderUI({
    rclipButton(id, "Copy to clipboard", text, icon("clipboard"))
  }))
}


shiny.huge.readCallTable <- function (fileName) {

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

  if (shiny.huge.isValidTable(ct)) {
    return(ct)
  } else {
    return(NULL)
  }

}

shiny.huge.isValidTable <- function (dt) {

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

shinyServer(function(input, output, session) {

  fullCallTable <- reactiveVal()

  sidebarFilteringReactives <- callModule(sidebarFiltering, "sidebarFiltering", fullCallTable, input$selectedColumns)

  geneTable <- sidebarFilteringReactives$geneTable
  filteredCallTable <- sidebarFilteringReactives$filteredCallTable
  comparedGenes <- reactiveVal()
  geneComparisonTable <- reactiveVal()

  dtInstance <- reactiveVal()

  ### CALL TABLE TAB

  observeEvent(input$callFile, {

    req(input$callFile)

    progress <- shiny::Progress$new()
    progress$set(message = "Uploading table", value = .5)

    ct <- shiny.huge.readCallTable(input$callFile$datapath)

    if(is.null(ct)) {
      shiny.huge.showErrorModal("Could not correctly read table. Make sure there are no trailing rows, all relevant columns are present and dots are used as decimal seperators.")
      progress$close()
    }

    req(ct)

    fullCallTable(ct)

    recognizedSymbolIndices <- shiny.huge.symbolToIndexMap[[ct$Symbol]]
    symbolsAreRecognized <- !is.na(recognizedSymbolIndices)
    recognizedSymbols <- unique(shiny.huge.geneTable$symbol[recognizedSymbolIndices[symbolsAreRecognized]])

    unexpressedSymbols <- recognizedSymbols[!recognizedSymbols %in% shiny.huge.gtexExpression$symbol]
    unrecognizedSymbols <- unique(ct$Symbol[!symbolsAreRecognized])

    unrecognizedSymbolsText <- paste0(unrecognizedSymbols, collapse = "\n")
    unexpressedSymbolsText <- paste0(unexpressedSymbols, collapse = "\n")

    symbolsInTotal <- length(unique(ct$Symbol))
    unrecognizedInTotal <- length(unrecognizedSymbols)
    unexpressedInTotal <- length(unexpressedSymbols)

    fractionUnrecognized <- round(unrecognizedInTotal / symbolsInTotal, digits = 2)
    fractionUnexpressed <- round(unexpressedInTotal / symbolsInTotal, digits = 2)

    output$unrecognizedSymbols <- renderText(unrecognizedSymbolsText)
    output$unexpressedSymbols <- renderText(unexpressedSymbolsText)

    output$totalUnrecognizedSymbols <- renderText(paste0(unrecognizedInTotal, " (", fractionUnrecognized, "%)"))
    output$totalUnexpressedSymbols <- renderText(paste0(unexpressedInTotal, " (", fractionUnexpressed, "%)"))

    output$unrecognizedToClipboard <- shiny.huge.copyToClipboardButton(unrecognizedSymbolsText, "unrecognizedbtn")
    output$unexpressedToClipboard <- shiny.huge.copyToClipboardButton(unexpressedSymbolsText, "unexpressedbtn")

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

  callModule(tableDownload, "filteredCallDownload", filteredCallTable, "filtered-calls-")

  ### GENE TABLE TAB

  output$geneTable <- DT::renderDataTable({
    req(geneTable())

    return(DT::datatable(geneTable(),
                         selection = "single",
                         style = "bootstrap",
                         class = DT:::DT2BSClass(c("hover", "stripe")))
    )
  })

  callModule(tableDownload, "geneDownload", geneTable, "genes-")

  ### ANNOTATION TABLE TAB

  output$annotationTable <- DT::renderDataTable(DT::datatable(shiny.huge.geneTable,
                         selection = "single",
                         style = "bootstrap",
                         class = DT:::DT2BSClass(c("compact", "hover", "stripe")),
                         options = list(columnDefs = list(list(
                           targets = seq_len(ncol(shiny.huge.geneTable)),
                           render = DT::JS(
                             "function(data, type, row, meta) {",
                             "if (data == null) return '';",
                             "return type === 'display' && data.length > 20 ?",
                             "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                             "}"
                           )))
                         )))

  callModule(tableDownload, "annotationDownload", function () shiny.huge.geneTable, "annotation-")
  callModule(tableDownload, "gtexDownload", function () shiny.huge.gtexExpression, "gtex-expression-")
  callModule(tableDownload, "hpaRnaDownload", function () shiny.huge.hpaRnaExpression, "hpa-rna-")
  callModule(tableDownload, "hpaProteinDownload", function () shiny.huge.hpaProteinExpession, "hpa-protein-")

  ### GENE COMPARISON TAB

  observeEvent(input$geneComparisonListUpload, {

    req(input$geneComparisonListUpload)

    comparedGenes(readLines(input$geneComparisonListUpload$datapath))

  })

  observe({

    req(comparedGenes())
    req(geneTable())

    if (input$geneComparisonUseFilters) {
      ourSymbols <- geneTable()$Symbol
    } else {
      ourSymbols <- unique(fullCallTable()$Symbol)
    }

    ourSymbols <- sort(ourSymbols)

    theirSymbols <- comparedGenes()

    ourIndices <- shiny.huge.symbolToIndexMap[[ourSymbols]]
    theirIndices <- shiny.huge.symbolToIndexMap[[theirSymbols]]

    ourIndexTable <- data.table(symbol = ourSymbols, index = ourIndices)
    theirIndexTable <- data.table(symbol = theirSymbols, index = theirIndices)

    ourNonNAIndices <- ourIndexTable[!is.na(index), index]
    theirNonNAIndices <- theirIndexTable[!is.na(index), index]

    onlyOurIndices <- setdiff(ourNonNAIndices, theirNonNAIndices)
    intersectingIndices <- intersect(ourNonNAIndices, theirNonNAIndices)
    onlyTheirIndices <- setdiff(theirNonNAIndices, ourNonNAIndices)

    onlyOurSymbols <- ourIndexTable[(index %in% onlyOurIndices | is.na(index)) & !symbol %in% theirIndexTable$symbol , symbol]
    intersectingSymbols <- ourIndexTable[index %in% intersectingIndices | symbol %in% theirIndexTable$symbol, symbol]
    onlyTheirSymbols <- theirIndexTable[(index %in% onlyTheirIndices | is.na(index)) & !symbol %in% ourIndexTable$symbol, symbol]

    onlyOurSymbolsText <- paste0(onlyOurSymbols, collapse = "\n")
    intersectingSymbolsText <- paste0(intersectingSymbols, collapse = "\n")
    onlyTheirSymbolsText <- paste0(onlyTheirSymbols, collapse = "\n")

    output$ourSymbols <- renderText(onlyOurSymbolsText)
    output$intersectingSymbols <- renderText(intersectingSymbolsText)
    output$theirSymbols <- renderText(onlyTheirSymbolsText)

    output$ourSymbolsToClipboard <- shiny.huge.copyToClipboardButton(onlyOurSymbolsText, "ourbtn")
    output$intersectingSymbolsToClipboard <- shiny.huge.copyToClipboardButton(intersectingSymbolsText, "intersectingbtn")
    output$theirSymbolsToClipboard <- shiny.huge.copyToClipboardButton(onlyTheirSymbolsText, "theirbtn")

    output$geneComparisonVenn <- renderPlot({

      vennInput <- list(
        our_symbols = c(shiny.huge.geneTable$symbol[ourNonNAIndices], ourIndexTable[is.na(index), symbol]),
        their_symbols = c(shiny.huge.geneTable$symbol[theirNonNAIndices], theirIndexTable[is.na(index), symbol])
      )

      venn <- venn.diagram(
        vennInput,
        filename = NULL,
        na = "remove",
        fontfamily = "sans",
        cat.fontfamily = "sans",
        fill = c("#3c8dbc", "orangered1")
      )

      grid.draw(venn)
    })
  })

  ### REACTOME TAB

  observeEvent(input$openReactomeButton, {

    req(geneTable())

    hgncSymbolIndices <- shiny.huge.symbolToIndexMap[[geneTable()$Symbol]]

    hgncGenes <- unique(shiny.huge.geneTable$symbol[hgncSymbolIndices[!is.na(hgncSymbolIndices)]])

    progress <- shiny::Progress$new()
    progress$set(message = paste("Opening reactome pathway overrepresenation link using", length(hgncGenes), "genes"), value = .5)

    apiUrl <- "https://reactome.org/AnalysisService/identifiers/projection"
    reactomeResponse <- POST(apiUrl, body = hgncGenes, add_headers("content-type" = "text/plain"), timeout(60))

    analysisToken <- content(reactomeResponse)$summary$token

    analysisUrl <- paste("https://reactome.org/PathwayBrowser/#DTAB=AN&ANALYSIS=", analysisToken, sep = "")

    tabOpeningJavascript <- paste("window.open('", analysisUrl, "', '_blank');", sep = "")

    output$reactomeOpeningScript <- renderUI(tags$script(HTML(tabOpeningJavascript)))

    progress$close()
  })

  ### DETAIL MODAL

  observeEvent(input$callTable_rows_selected,
               shiny.huge.geneExpressionModal(
                 filteredCallTable()[input$callTable_rows_selected]$Symbol,
                 filteredCallTable,
                 input,
                 output)
  )

  observeEvent(input$geneTable_rows_selected,
               shiny.huge.geneExpressionModal(
                 geneTable()[input$geneTable_rows_selected, Symbol],
                 filteredCallTable,
                 input,
                 output)
  )

  observeEvent(input$modalOkBtn, {
    removeModal()
  })


})
