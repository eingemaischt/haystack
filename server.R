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

    errorMessage <- "No gtex information available"

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
      ylim(c(0, NA)) +
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

  actualCharacterColumnTypes <- sapply(dt[,..characterColumns], class)
  actualNumericColumnTypes <- sapply(dt[,..numericColumns], class)

  return(
    all(expectedColumns %in% colnames(dt)) &&
    all(actualCharacterColumnTypes == "character") &&
    all(actualNumericColumnTypes %in% c("numeric", "integer"))
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
        tabPanel("GTEx tissues (TPM)", plotOutput("modalGTExExpression", height = "640px")),
        tabPanel("GTex tissues (TPM scaled)", plotOutput("modalGTExScaledExpression", height = "640px")),
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

    expression <- shiny.huge.gtexExpression[
      symbol %in% matchingGene$symbol &
      tpm >= input$minRawTPM
    ]

    rawValues <- expression[,list(tissue = tissue, value = tpm)]
    scaledValues <- expression[,list(tissue = tissue, value = tpm_scaled)]

    callsInGene <- callTableReactiveVal()[Symbol == selectedSymbol]
    # Patients may be sampled twice, so we need to 'pre-deduplicate' before the
    # actual deduplication
    callsPerPatientInGene <- callsInGene[,
      list("Read depth" = mean(`Read depth`), "Variant depth" = mean(`Variant depth`), "Samples" = .N),
      by = list(HGVSc, AlternativePatNr, Chr, Position, Consequence, `AF Popmax`)
    ]
    uniqueMutations <- callsPerPatientInGene[,
      list("Mean read depth" = mean(`Read depth`), "Mean variant depth" = mean(`Variant depth`), "Patients" = .N, "Samples" = sum(Samples)),
      by = list(HGVSc, Chr, Position, Consequence, `AF Popmax`)
    ]

    output$modalGTExExpression <- shiny.huge.modalExpressionPlot(rawValues, input$expressions, paste(selectedSymbol, "GTEx data (raw TPM)", sep = ": "))
    output$modalGTExScaledExpression <- shiny.huge.modalExpressionPlot(scaledValues, input$expressions, paste(selectedSymbol, "GTEx data (scaled TPM)", sep = ": "))
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

shiny.huge.resetFilters <- function (session, callTable) {

  numberOfPatients <- length(unique(callTable$AlternativePatNr))
  consequences <- unique(unlist(strsplit(callTable$Consequence, ",")))

  updateCheckboxGroupInput(session, "selectedColumns", choices = colnames(callTable), selected = colnames(callTable))
  updateSliderInput(session, "patientNumber", value = c(0, numberOfPatients), max = numberOfPatients)
  updateSliderInput(session, "minReadDepth", value = 0, max = max(callTable$`Read depth`))
  updateSliderInput(session, "minVariantDepth", value = 0, max = max(callTable$`Variant depth`, na.rm = TRUE))
  updateSliderInput(session, "scaledTPM", value = c(0,1))
  updateCheckboxGroupInput(session, "genotypes", selected = c("unknown", "hom_ref", "het", "hom_alt"))
  updateCheckboxInput(session, "onlyCompoundHeterozygosity", value = FALSE)
  updateNumericInput(session, "maxAFPopmax", value = 100)
  updateNumericInput(session, "minRawTPM", value = 0)
  updateSelectizeInput(session, "expressions", selected = NULL, choices = unique(shiny.huge.gtexExpression$tissue))
  updateSelectizeInput(session, "consequences", selected = NULL, choices = consequences)

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

    recognizedSymbols <- shiny.huge.geneTable$symbol[shiny.huge.symbolToIndexMap[[ct$Symbol]]]

    unrecognizedCalls <- ct[!recognizedSymbols %in% shiny.huge.gtexExpression$symbol]
    unrecognizedSymbols <- unique(unrecognizedCalls$Symbol)

    output$unrecognizedSymbols <- renderText(paste0(unrecognizedSymbols, collapse = "\n"))

    progress$close()
  })

  observe({

    maxPopMax <- input$maxAFPopmax
    minRawTPM <- input$minRawTPM

    shiny.huge.handleErrorNotification(maxPopMax, "AF Popmax ", "popmaxErrorNotification")
    shiny.huge.handleErrorNotification(minRawTPM, "minimum TPM", "minRawTPMErrorNotification")

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

  ### GENE TABLE TAB

  output$geneTable <- DT::renderDataTable({
    req(geneTable())

    return(DT::datatable(geneTable(),
                         selection = "single",
                         style = "bootstrap",
                         class = DT:::DT2BSClass(c("hover", "stripe")))
    )
  })

  output$geneDownload <- downloadHandler(
    filename = "genes.csv",
    content = function (con) write.csv(geneTable(), file = con, row.names = FALSE)
  )

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

    output$ourSymbols <- renderText(paste0(onlyOurSymbols, collapse = "\n"))
    output$intersectingSymbols <- renderText(paste0(intersectingSymbols, collapse = "\n"))
    output$theirSymbols <- renderText(paste0(onlyTheirSymbols, collapse = "\n"))

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

    hgncGenes <- unique(shiny.huge.geneTable$symbol[hgncSymbolIndices])

    progress <- shiny::Progress$new()
    progress$set(message = paste("Opening reactome pathway overrepresenation link using", length(hgncGenes), "genes"), value = .5)

    apiUrl <- "https://reactome.org/AnalysisService/identifiers/projection"
    reactomeResponse <- POST(apiUrl, body = hgncGenes, add_headers("content-type" = "text/plain"))

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
    req(input$patientNumber)
    req(input$minReadDepth)
    req(input$minVariantDepth)
    req(input$maxAFPopmax)
    req(input$minRawTPM)
    req(input$scaledTPM)

    ct <- fullCallTable()

    patientSymbolStrings <- paste0(ct$AlternativePatNr, ct$Symbol)
    patientSymbolDuplicates <- duplicated(patientSymbolStrings) | duplicated(patientSymbolStrings, fromLast = TRUE)

    ct <- ct[
      (!input$onlyCompoundHeterozygosity | patientSymbolDuplicates) &
        (
          ("hom_ref" %in% input$genotypes & Genotype == "0/0") |
          ("het" %in% input$genotypes & (Genotype == "0/1" | Genotype == "1/0")) |
          ("hom_alt" %in% input$genotypes & Genotype == "1/1") |
          ("unknown" %in% input$genotypes & grepl(".", Genotype, fixed = TRUE))
        ) &
        `Read depth` >= input$minReadDepth &
        (`Variant depth` >= input$minVariantDepth | is.na(`Variant depth`)) &
        (is.na(`AF Popmax`) | input$maxAFPopmax >= `AF Popmax`) &
        (is.null(input$consequences) | grepl(paste0(input$consequences, collapse = "|"), Consequence)),
      input$selectedColumns, with = FALSE
      ]

    # count mutations per gene per patient by multiple aggregations
    gt <- ct[, list(patients = .N), by = .(Symbol, AlternativePatNr)]
    gt <- gt[, list(patients = .N, compound_het_patients = sum(patients > 1)), by = .(Symbol)]
    gt <- gt[order(patients, decreasing = TRUE)]

    matchingGeneIndices <- shiny.huge.symbolToIndexMap[[gt$Symbol]]
    matchingGeneNames <- shiny.huge.geneTable$symbol[matchingGeneIndices]

    gt$name <- shiny.huge.geneTable$name[matchingGeneIndices]
    gt$locus_type <- shiny.huge.geneTable$locus_type[matchingGeneIndices]
    gt$family <- shiny.huge.geneTable$gene_family[matchingGeneIndices]
    gt$description <- shiny.huge.geneTable$description[matchingGeneIndices]

    expressionFilteredGenes <- shiny.huge.gtexExpression[
      tpm >= input$minRawTPM &
      tpm_scaled >= input$scaledTPM[1] & tpm_scaled <= input$scaledTPM[2] &
      tissue %in% input$expressions
    ]

    gt <- gt[
      input$patientNumber[1] <= patients &
      input$patientNumber[2] >= patients &
      (is.null(input$expressions) | matchingGeneNames %in% expressionFilteredGenes$symbol)
      ]
    ct <- ct[Symbol %in% gt$Symbol]

    geneTable(gt)
    filteredCallTable(ct)

  })

})
