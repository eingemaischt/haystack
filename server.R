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

shiny.huge.handleErrorNotification <- function (value, filterName, notificationId) {

  if (is.na(value)) {
    showNotification(paste0("Invalid number in ", filterName, " filter"), closeButton = FALSE, type = "error", duration = NULL, id = notificationId)
  } else {
    removeNotification(notificationId)
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

shiny.huge.readCallTable <- function (fileName) {

  ct <- fread(fileName)

  if (is.numeric(ct$Chr)) {
    ct$Chr <- as.character(ct$Chr)
  }

  if (shiny.huge.isValidTable(ct)) {
    return(ct)
  } else {
    return(NULL)
  }

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
        tabPanel("GTEx (TPM)", plotOutput("modalGTExExpression", height = "640px")),
        tabPanel("GTex (TPM scaled)", plotOutput("modalGTExScaledExpression", height = "640px")),
        tabPanel("HPA RNA (TPM)", plotOutput("modalHpaRnaExpression", height = "640px")),
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
      list("Mean read depth" = mean(`Read depth`), "Mean variant depth" = mean(`Variant depth`), "Samples" = length(unique(Sample))),
      by = list(HGVSc, Chr, Position, Consequence, `AF Popmax`)
    ]

    matchingGtexExpression <- shiny.huge.gtexExpression[symbol == matchingGene$symbol]
    matchingHpaRnaExpression <- shiny.huge.hpaRnaExpression[symbol == matchingGene$symbol]
    matchingHpaProteinExpression <- shiny.huge.hpaProteinExpession[symbol == matchingGene$symbol]

    rawGtexValues <- matchingGtexExpression[,list(tissue = tissue, value = tpm)]
    scaledGtexValues <- matchingGtexExpression[,list(tissue = tissue, value = tpm_scaled)]

    rawHpaRnaValues <- matchingHpaRnaExpression[,list(tissue = tissue, value = tpm)]
    scaledHpaRnaValues <- matchingHpaRnaExpression[,list(tissue = tissue, value = tpm_scaled)]

    hpaProteinValues <- matchingHpaProteinExpression[,list(value = max(level)),by = tissue]

    output$modalGTExExpression <- shiny.huge.modalExpressionPlot(rawGtexValues, input$expressions, paste(selectedSymbol, "GTEx data (raw TPM)", sep = ": "))
    output$modalGTExScaledExpression <- shiny.huge.modalExpressionPlot(scaledGtexValues, input$expressions, paste(selectedSymbol, "GTEx data (scaled TPM)", sep = ": "))
    output$modalHpaRnaExpression <- shiny.huge.modalExpressionPlot(rawHpaRnaValues, input$expressions, paste(selectedSymbol, "HPA RNA data (raw TPM)", sep = ":"))
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

shiny.huge.handleTableDownload <- function (tableReactiveValue, filePrefix) {

  return(downloadHandler(
    filename = function() {
      paste0(filePrefix, Sys.time(), ".csv")
    },
    content = function (con) {

      if(is.null(tableReactiveValue())) {

        writtenData <- data.table("NO DATA AVAILABLE")

      } else {

        writtenData <- tableReactiveValue()

      }

      fwrite(writtenData, con)
    }
  ))

}

shiny.huge.resetFilters <- function (session, callTable) {

  numberOfSamples <- length(unique(callTable$Sample))
  consequences <- unique(unlist(strsplit(callTable$Consequence, ",")))

  studies <- unique(callTable$Studie)

  updateCheckboxGroupInput(session, "selectedColumns", choices = colnames(callTable), selected = colnames(callTable))
  updateSliderInput(session, "sampleNumber", value = c(0, numberOfSamples), max = numberOfSamples)
  updateSliderInput(session, "minReadDepth", value = 0)
  updateSliderInput(session, "minVariantDepth", value = 0)
  updateSliderInput(session, "readVariantFrequency", value = c(0,1))
  updateSliderInput(session, "scaledTPM", value = c(0,1))
  updateSelectInput(session, "proteinLevel", selected = "Any")
  updateCheckboxGroupInput(session, "genotypes", selected = c("unknown", "het", "hom_alt"))
  updateCheckboxInput(session, "onlyCompoundHeterozygosity", value = FALSE)
  updateNumericInput(session, "maxAFPopmax", value = 100)
  updateSelectizeInput(session, "expressions", selected = NULL, choices = unique(c(shiny.huge.gtexExpression$tissue, shiny.huge.hpaRnaExpression$tissue, shiny.huge.hpaProteinExpession$tissue)))
  updateSelectizeInput(session, "consequences", selected = NULL, choices = consequences)
  updateSelectizeInput(session, "studies", selected = NULL, choices = studies)
  updateSelectizeInput(session, "chromosomes", selected = NULL, choices = unique(callTable$Chr))

}

shiny.huge.copyToClipboardButton <- function (text, id) {
  return(renderUI({
    rclipButton(id, "Copy to clipboard", text, icon("clipboard"))
  }))
}

shinyServer(function(input, output, session) {

  fullCallTable <- reactiveVal()
  filteredCallTable <- reactiveVal()

  geneTable <- reactiveVal()
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

    shiny.huge.resetFilters(session, ct)

    fullCallTable(ct)

    recognizedSymbolIndices <- shiny.huge.symbolToIndexMap[[ct$Symbol]]
    symbolsAreRecognized <- !is.na(recognizedSymbolIndices)
    recognizedSymbols <- unique(shiny.huge.geneTable$symbol[recognizedSymbolIndices[symbolsAreRecognized]])

    unexpressedSymbols <- recognizedSymbols[!recognizedSymbols %in% unique(c(shiny.huge.gtexExpression$symbol, shiny.huge.hpaRnaExpression$symbol, shiny.huge.hpaProteinExpession$symbol))]
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

  observe({

    maxPopMax <- input$maxAFPopmax

    shiny.huge.handleErrorNotification(maxPopMax, "AF Popmax ", "popmaxErrorNotification")

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

  output$filteredCallTableDownload <- shiny.huge.handleTableDownload(filteredCallTable, "filtered-calls-")

  ### GENE TABLE TAB

  output$geneTable <- DT::renderDataTable({
    req(geneTable())

    return(DT::datatable(geneTable(),
                         selection = "single",
                         style = "bootstrap",
                         class = DT:::DT2BSClass(c("hover", "stripe")))
    )
  })

  output$geneDownload <- shiny.huge.handleTableDownload(geneTable, "genes-")

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

  output$annotationDownload <- shiny.huge.handleTableDownload(function () shiny.huge.geneTable, "annotation-")
  output$gtexDownload <- shiny.huge.handleTableDownload(function () shiny.huge.gtexExpression, "gtex-expression-")
  output$hpaRnaDownload <- shiny.huge.handleTableDownload(function () shiny.huge.hpaRnaExpression, "hpa-rna-")
  output$hpaProteinDownload <- shiny.huge.handleTableDownload(function () shiny.huge.hpaProteinExpession, "hpa-protein-")

  ### GENE COMPARISON TAB

  observeEvent(input$geneComparisonListUpload, {

    req(input$geneComparisonListUpload)

    comparedGenes(readLines(input$geneComparisonListUpload$datapath))

  })

  observe({

    req(comparedGenes())
    req(geneTable())

    ourSymbols <- geneTable()$Symbol
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

      venn <- venn.diagram(list(
        our_symbols = c(shiny.huge.geneTable$symbol[ourNonNAIndices], ourIndexTable[is.na(index), symbol]),
        their_symbols = c(shiny.huge.geneTable$symbol[theirNonNAIndices], theirIndexTable[is.na(index), symbol])
      ), filename = NULL, na = "remove")

      grid.draw(venn)
    })
  })

  ### SIDEBAR

  output$numTotalRows    <- renderText({ paste("Total calls in table: ",    nrow(fullCallTable())) })
  output$numFilteredRows <- renderText({ paste("Filtered calls in table: ", nrow(filteredCallTable())) })

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

  ### FILTERING

  observeEvent(input$filterReset, {

    req(fullCallTable())

    shiny.huge.resetFilters(session, fullCallTable())

  })

  observe({

    req(fullCallTable())
    req(input$selectedColumns)
    req(input$sampleNumber)
    req(input$minReadDepth)
    req(input$minVariantDepth)
    req(input$maxAFPopmax)
    req(input$scaledTPM)
    req(input$proteinLevel)

    ct <- fullCallTable()

    variantFrequency <- ct$`Variant depth` / ct$`Read depth`

    ct <- ct[
      (
        ("het" %in% input$genotypes & (Genotype == "0/1" | Genotype == "1/0")) |
        ("hom_alt" %in% input$genotypes & Genotype == "1/1") |
        ("unknown" %in% input$genotypes & (grepl(".", Genotype, fixed = TRUE) | Genotype == "0/0"))
      ) &
      `Read depth` >= input$minReadDepth &
      (`Variant depth` >= input$minVariantDepth | is.na(`Variant depth`)) &
      ((input$readVariantFrequency[1] <= variantFrequency & variantFrequency <= input$readVariantFrequency[2]) | `Read depth` == 0) &
      (is.na(`AF Popmax`) | input$maxAFPopmax >= `AF Popmax`) &
      (is.null(input$consequences) | grepl(paste0(input$consequences, collapse = "|"), Consequence)) &
      (is.null(input$studies) | grepl(paste0(input$studies, collapse = "|"), Studie)) &
      (is.null(input$chromosomes) | Chr %in% input$chromosomes),
      input$selectedColumns, with = FALSE
    ]

    # Filter for compound heterozygosity only after other variant filters have been applied
    # Otherwise, there will be variants that are not compound based on the filters
    sampleSymbolStrings <- paste0(ct$Sample, ct$Symbol)
    sampleSymbolDuplicates <- duplicated(sampleSymbolStrings) | duplicated(sampleSymbolStrings, fromLast = TRUE)

    ct <- ct[!input$onlyCompoundHeterozygosity | sampleSymbolDuplicates]

    # count mutations per gene per sample by multiple aggregations
    gt <- ct[, list(samples = .N), by = .(Symbol, Sample)]
    gt <- gt[, list(samples = .N, compound_het_samples = sum(samples > 1)), by = .(Symbol)]
    gt <- gt[order(samples, decreasing = TRUE)]

    matchingGeneIndices <- shiny.huge.symbolToIndexMap[[gt$Symbol]]
    matchingGeneNames <- shiny.huge.geneTable$symbol[matchingGeneIndices]

    gt$name <- shiny.huge.geneTable$name[matchingGeneIndices]
    gt$locus_type <- shiny.huge.geneTable$locus_type[matchingGeneIndices]
    gt$family <- shiny.huge.geneTable$gene_family[matchingGeneIndices]
    gt$description <- shiny.huge.geneTable$description[matchingGeneIndices]

    gtexFilteredGenes <- shiny.huge.gtexExpression[
      tpm_scaled >= input$scaledTPM[1] & tpm_scaled <= input$scaledTPM[2] &
      tissue %in% input$expressions
    ]

    hpaRnaFilteredGenes <- shiny.huge.hpaRnaExpression[
      tpm_scaled >= input$scaledTPM[1] & tpm_scaled <= input$scaledTPM[2] &
      tissue %in% input$expressions
    ]

    hpaProteinFilteredGenes <- shiny.huge.hpaProteinExpession[
      input$proteinLevel <= level &
      tissue %in% input$expressions
    ]

    gt <- gt[
      input$sampleNumber[1] <= samples &
      input$sampleNumber[2] >= samples &
      (is.null(input$expressions) | matchingGeneNames %in% c(unique(gtexFilteredGenes$symbol), unique(hpaRnaFilteredGenes$symbol))) &
      (is.null(input$expressions) | input$proteinLevel == "Any" | matchingGeneNames %in% hpaProteinFilteredGenes$symbol)
      ]
    ct <- ct[Symbol %in% gt$Symbol]

    geneTable(gt)
    filteredCallTable(ct)

  })

})
