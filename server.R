options(shiny.maxRequestSize=100*1024^2)

shiny.huge.detailDivElement <- function (label, value) {

  return(tags$div(
    tags$b(label),
    tags$em(value)
  ))

}

shiny.huge.modalExpressionPlot <- function (expressions, expressionFilter, title) {

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

shiny.huge.geneModalExpression <- function (callTableReactiveVal, row_selection, input, output) {

  return({

    req(callTableReactiveVal())

    selectedRow <- callTableReactiveVal()[row_selection]

    matchingGene <- shiny.huge.geneTable[
      symbol == selectedRow$Symbol |
        grepl(selectedRow$Symbol, prev_symbol) |
        grepl(selectedRow$Symbol, alias_symbol)
      ]

    showModal(modalDialog(
      shiny.huge.detailDivElement("Location:", matchingGene$location),
      shiny.huge.detailDivElement("Gene (found in table):", selectedRow$Symbol),
      shiny.huge.detailDivElement("Gene (current HGNC symbol):", matchingGene$symbol),
      shiny.huge.detailDivElement("Gene full name:", matchingGene$name),
      shiny.huge.detailDivElement("Gene family:", matchingGene$gene_family),
      shiny.huge.detailDivElement("Gene description:", matchingGene$description),
      shiny.huge.detailDivElement("Location type:", matchingGene$location_type),
      shiny.huge.detailDivElement("Ensembl ID:", matchingGene$ensembl_gene_id),
      tags$hr(),
      tags$b("Expression:"),
      tabsetPanel(
        tabPanel("GTEx tissues", plotOutput("modalGTExExpression", height = "640px"))
      ),
      title = "Details",
      footer = actionButton("modalOkBtn", label = "OK", icon = icon("ok")),
      size = "l",
      easyClose = TRUE
    ))

    expression <- shiny.huge.queryExpressionAnnotations(selectedRow$Symbol)[[1]]
    # workaround for 'max' not working properly for factors in data.table
    # see: https://github.com/Rdatatable/data.table/issues/1947
    maxFixFn <- max
    expression <- lapply(expression, function (expr) {

      if (nrow(expr) == 0) return(expr)

      expr[,list(value = maxFixFn(value)), by = "tissue"]
    })

    #output$modalNormalExpression <- shiny.huge.modalExpressionPlot(expression$normal, input$expressions, paste(selectedRow$RSID, "HPA immunohistochemistry data", sep = ": "))
    #output$modalRnaExpression <- shiny.huge.modalExpressionPlot(expression$rna, input$expressions, paste(selectedRow$RSID, "HPA RNA-seq data", sep = ": "))
    output$modalGTExExpression <- shiny.huge.modalExpressionPlot(expression$gtex, input$expressions, paste(selectedRow$Symbol, "GTEx data", sep = ": "))
  })

}

shinyServer(function(input, output, session) {

  fullCallTable <- reactiveVal()
  filteredCallTable <- reactiveVal()

  geneTable <- reactiveVal()

  dtInstance <- reactiveVal()

  ### CALL TABLE TAB

  observeEvent(input$callFile, {

    req(input$callFile)

    progress <- shiny::Progress$new()
    progress$set(message = "Uploading table", value = .5)

    ct <- fread(input$callFile$datapath)

    req(ct)

    updateCheckboxGroupInput(session, "selectedColumns", choices = colnames(ct), selected = colnames(ct))
    updateNumericInput(session, "minSamplePercentage", value = 0)
    updateSelectizeInput(session, "expressions", selected = NULL)

    fullCallTable(ct)

    progress$close()
  })

  observe({

    req(fullCallTable())
    req(input$selectedColumns)

    ct <- fullCallTable()

    ct <- ct[,input$selectedColumns, with = FALSE]

    # count mutations per gene per sample by multiple aggregations
    gt <- ct[, .N, by = .(Symbol, Sample)]
    gt <- gt[, .N, by = .(Symbol)]
    gt <- gt[order(N, decreasing = TRUE)]

    matchingGeneIndices <- match(gt$Symbol, shiny.huge.geneTable$symbol)

    gt$name <- shiny.huge.geneTable$name[matchingGeneIndices]
    gt$locus_type <- shiny.huge.geneTable$locus_type[matchingGeneIndices]
    gt$family <- shiny.huge.geneTable$gene_family[matchingGeneIndices]

    geneTable(gt)
    filteredCallTable(ct)
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

  output$geneTable <- DT::renderDataTable({
    req(geneTable())

    return(DT::datatable(geneTable(),
                         selection = "single",
                         style = "bootstrap",
                         class = DT:::DT2BSClass(c("hover", "stripe")))
    )
  })

  ### SIDEBAR

  output$numTotalRows    <- renderText({ paste("Total calls in table: ",    nrow(fullCallTable())) })
  output$numFilteredRows <- renderText({ paste("Filtered calls in table: ", nrow(filteredCallTable())) })

  ### REACTOME TAB

  observeEvent(input$openReactomeButton, {

    req(filteredCallTable())
    req(filteredCallTable()$symbol)

    hgncGenes <- unique(filteredCallTable()$symbol)
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

  observeEvent(input$callTable_rows_selected,
               shiny.huge.geneModalExpression(
                 filteredCallTable,
                 input$callTable_rows_selected,
                 input,
                 output)
  )

  observeEvent(input$geneTable_rows_selected,
               shiny.huge.geneModalExpression(
                 geneTable,
                 input$geneTable_rows_selected,
                 input,
                 output)
  )

  observeEvent(input$modalOkBtn, {
    removeModal()
  })
})
