options(shiny.maxRequestSize=100*1024^2)

shiny.p7.histogramExpression <- function (reactiveTableValue, colName, statType = "bin", numBins = 50) {
  return({
    req(reactiveTableValue())
    req(colName %in% colnames(reactiveTableValue()))

    ggplot(reactiveTableValue(), aes_string(x = colName)) +
      geom_histogram(bins = numBins, stat = statType, fill = "steelblue", colour = "steelblue4") +
      ggtitle(paste(colName, "histogram"))
  })
}

shiny.p7.detailDivElement <- function (label, value) {

  return(tags$div(
    tags$b(label),
    tags$em(value)
  ))

}

shiny.p7.modalExpressionPlot <- function (expressions, expressionFilter, title) {

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

shiny.p7.downloadHandler <- function (tableReactiveVal, fileName) {

  return(downloadHandler(
    filename = fileName,
    content = function (file) {
      write.table(tableReactiveVal(), file = file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
  ))

}

shiny.p7.applyCandidateModalExpression <- function (candidateValue, rowsSelected, snpTable, fullTable, proxy) {

  selectedRow <- snpTable()[rowsSelected]

  fullDt <- fullTable()
  fullDt[ RSID == selectedRow$RSID, IS_CANDIDATE := candidateValue ]

  snpDt <- snpTable()
  snpDt[ RSID == selectedRow$RSID, IS_CANDIDATE := candidateValue ]

  replaceData(proxy, snpDt, resetPaging = FALSE, clearSelection = TRUE)

  removeModal()
}

shinyServer(function(input, output, session) {

  fullTable <- reactiveVal()
  snpTable <- reactiveVal()
  dtInstance <- reactiveVal()

  ### SNP TABLE TAB

  observeEvent(input$snpFile, {

    req(input$snpFile)

    progress <- shiny::Progress$new()
    progress$set(message = "Querying SNP metadata from databases", value = .5)

    dt <- tryCatch({
      shiny.p7.annotateSNPFile(input$snpFile$datapath)
    }, error = function (e) {
      showModal(modalDialog(title = "Biomart not reachable", tags$p("The ensembl server is down (again...), so please try again later.")))
      progress$close()
      return(NULL)
    })

    req(dt)

    updateCheckboxGroupInput(session, "selectedColumns", choices = colnames(dt), selected = colnames(dt))

    possibleExpressions <- unique(unlist(strsplit(dt$EXPRESSIONS, ",")))
    updateSelectizeInput(session, "expressions", choices = unique(possibleExpressions), selected = NULL)

    fullTable(dt)

    progress$close()
  })

  observe({

    req(fullTable())
    req(input$selectedColumns)
    req(input$maxPValue)
    req(input$minAbsBeta)
    req(input$minNMiss)
    req(input$minR2)
    req(input$emptyGenesRadio)

    dt <- fullTable()
    dt <- dt[
      P <=
        input$maxPValue &
        abs(BETA) >= input$minAbsBeta &
        NMISS >= input$minNMiss &
        R2 >= input$minR2 &
        (input$emptyGenesRadio == "TRUE" | nchar(ENSEMBL_GENES) > 0) &
        (is.null(input$expressions) | (input$emptyGenesRadio == "TRUE" & nchar(ENSEMBL_GENES) == 0) | seq_len(nrow(dt)) %in% grep(paste(input$expressions, collapse = "|"), EXPRESSIONS)) &
        (is.null(input$candidateFilter) | seq_len(nrow(dt)) %in% grep(paste(input$candidateFilter, collapse = "|"), IS_CANDIDATE)),
      input$selectedColumns, with = FALSE]

    snpTable(dt)
  })

  output$snpTable <- DT::renderDataTable({

    req(snpTable())
    return(DT::datatable(snpTable(),
                         selection = "single",
                         style = "bootstrap",
                         class = DT:::DT2BSClass(c("compact", "hover", "stripe")),
                         options = list(columnDefs = list(list(
                          targets = seq_len(ncol(snpTable())),
                          render = DT::JS(
                            "function(data, type, row, meta) {",
                            "if (data == null) return '';",
                            "return type === 'display' && data.length > 20 ?",
                            "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                            "}"
                          )))
                        )))
  })


  output$downloadAllData <- shiny.p7.downloadHandler(fullTable, "whole-table.csv")
  output$downloadFilteredData <- shiny.p7.downloadHandler(snpTable, "filtered-table.csv")

  ### SIDEBAR

  output$numTotalRows    <- renderText({ paste("Total SNPs in table: ",    nrow(fullTable())) })
  output$numFilteredRows <- renderText({ paste("Filtered SNPs in table: ", nrow(snpTable())) })

  ### HISTOGRAM TAB

  output$betaHistPlot <- renderPlot(shiny.p7.histogramExpression(snpTable, "BETA"),  res = 90)
  output$pHistPlot    <- renderPlot(shiny.p7.histogramExpression(snpTable, "P"),     res = 90)
  output$r2Plot       <- renderPlot(shiny.p7.histogramExpression(snpTable, "R2"),    res = 90)
  output$nmissPlot    <- renderPlot(shiny.p7.histogramExpression(snpTable, "NMISS"), res = 90)

  ### REACTOME TAB

  observeEvent(input$openReactomeButton, {

    req(snpTable())
    req(snpTable()$HGNC_NAMES)

    hgncGenes <- strsplit(snpTable()$HGNC_NAMES, split = ",")
    hgncGenes <- unlist(hgncGenes)
    hgncGenes <- unique(hgncGenes)
    hgncGenes <- hgncGenes[hgncGenes != "NULL"]

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

  observeEvent(input$snpTable_rows_selected, {

    req(snpTable())

    selectedRow <- snpTable()[input$snpTable_rows_selected]
    isCandidate = selectedRow$IS_CANDIDATE

    showModal(modalDialog(
      shiny.p7.detailDivElement("SNP identifier:", selectedRow$RSID),
      shiny.p7.detailDivElement("Chromosome:", selectedRow$CHR),
      shiny.p7.detailDivElement("Genes:", selectedRow$HGNC_NAMES),
      shiny.p7.detailDivElement("Gene descriptions:", selectedRow$GENE_DESCRIPTION),
      shiny.p7.detailDivElement("Ensembl IDs:", selectedRow$ENSEMBL_GENES),
      actionButton("modalApproveCandidateBtn",
                   label = "Approve candidate",
                   icon = icon("thumbs-up"),
                   class = "btn btn-sm btn-success"),
      actionButton("modalDiscardCandidateBtn",
                   label = "Discard candidate",
                   icon = icon("thumbs-down"),
                   class = "btn btn-sm btn-danger"),
      actionButton("modalResetCandidateBtn",
                   label = "Reset candidate",
                   icon = icon("asterisk"),
                   class = "btn btn-sm btn-default"),
      tags$hr(),
      tags$b("Expression:"),
      tabsetPanel(
        tabPanel("HPA normal tissues", plotOutput("modalNormalExpression", height = "640px")),
        tabPanel("HPA RNA tissues", plotOutput("modalRnaExpression", height = "640px")),
        tabPanel("GTEx tissues", plotOutput("modalGTExExpression", height = "640px"))
      ),
      title = "Details",
      footer = actionButton("modalOkBtn", label = "OK", icon = icon("ok")),
      size = "l",
      easyClose = TRUE
    ))

    expression <- shiny.p7.queryExpressionAnnotations(selectedRow$ENSEMBL_GENES)[[1]]
    # workaround for 'max' not working properly for factors in data.table
    # see: https://github.com/Rdatatable/data.table/issues/1947
    maxFixFn <- max
    expression <- lapply(expression, function (expr) {

      if (nrow(expr) == 0) return(expr)

      expr[,list(value = maxFixFn(value)), by = "tissue"]
    })

    output$modalNormalExpression <- shiny.p7.modalExpressionPlot(expression$normal, input$expressions, paste(selectedRow$RSID, "HPA immunohistochemistry data", sep = ": "))
    output$modalRnaExpression <- shiny.p7.modalExpressionPlot(expression$rna, input$expressions, paste(selectedRow$RSID, "HPA RNA-seq data", sep = ": "))
    output$modalGTExExpression <- shiny.p7.modalExpressionPlot(expression$gtex, input$expressions, paste(selectedRow$RSID, "GTEx data", sep = ": "))
  })

  proxy <- dataTableProxy("snpTable")
  observeEvent(input$modalApproveCandidateBtn, {

    shiny.p7.applyCandidateModalExpression("approved", input$snpTable_rows_selected, snpTable, fullTable, proxy)

  })

  observeEvent(input$modalDiscardCandidateBtn, {

    shiny.p7.applyCandidateModalExpression("discarded", input$snpTable_rows_selected, snpTable, fullTable, proxy)

  })

  observeEvent(input$modalResetCandidateBtn, {

    shiny.p7.applyCandidateModalExpression("none", input$snpTable_rows_selected, snpTable, fullTable, proxy)

  })

  observeEvent(input$modalOkBtn, {
    removeModal()
  })
})
