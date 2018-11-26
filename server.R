shinyServer(function(input, output, session) {

  callTableReactives <- callModule(callTable, "callTableTab")

  fullCallTable <- callTableReactives$fullCallTable
  filteredCallTable <- callTableReactives$filteredCallTable

  geneTable <- callModule(sidebarFiltering, "sidebarFiltering", fullCallTable, filteredCallTable, callTableReactives$selectedColumns)

  callModule(geneTableTab, "geneTableTab", geneTable, filteredCallTable)
  callModule(annotationsTab, "annotationsTab")

  comparedGenes <- reactiveVal()
  geneComparisonTable <- reactiveVal()



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

    output$ourSymbolsToClipboard <- copyToClipboardButton(onlyOurSymbolsText, "ourbtn")
    output$intersectingSymbolsToClipboard <- copyToClipboardButton(intersectingSymbolsText, "intersectingbtn")
    output$theirSymbolsToClipboard <- copyToClipboardButton(onlyTheirSymbolsText, "theirbtn")

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


})
