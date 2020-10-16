tabs.geneComparison.module <- function (input, output, session, fullCallTable, geneTable) {

  comparedGenes <- reactiveVal()
  geneComparisonTable <- reactiveVal()

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

    ourIndices <- annotation.symbolToIndexMap[ourSymbols]
    theirIndices <- annotation.symbolToIndexMap[theirSymbols]

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

    output$ourSymbolsToClipboard <-util.copyToClipboardButton(onlyOurSymbolsText, "ourbtn")
    output$intersectingSymbolsToClipboard <-util.copyToClipboardButton(intersectingSymbolsText, "intersectingbtn")
    output$theirSymbolsToClipboard <-util.copyToClipboardButton(onlyTheirSymbolsText, "theirbtn")

    output$geneComparisonVenn <- renderPlot({

      vennInput <- list(
        our_symbols = c(annotation.geneTable$symbol[ourNonNAIndices], ourIndexTable[is.na(index), symbol]),
        their_symbols = c(annotation.geneTable$symbol[theirNonNAIndices], theirIndexTable[is.na(index), symbol])
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

}
