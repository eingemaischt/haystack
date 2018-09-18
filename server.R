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

shinyServer(function(input, output, session) {

  fullTable <- reactiveVal()
  filteredTable <- reactiveVal()
  dtInstance <- reactiveVal()

  ### CALL TABLE TAB

  observeEvent(input$callFile, {

    req(input$callFile)

    progress <- shiny::Progress$new()
    progress$set(message = "Uploading table", value = .5)

    dt <- fread(input$callFile$datapath)

    req(dt)

    updateCheckboxGroupInput(session, "selectedColumns", choices = colnames(dt), selected = colnames(dt))

    updateSelectizeInput(session, "expressions", choices = unique(shiny.huge.gtexExpression$tissue), selected = NULL)

    fullTable(dt)

    progress$close()
  })

  observe({

    req(fullTable())
    req(input$selectedColumns)

    dt <- fullTable()
    dt <- dt[,input$selectedColumns, with = FALSE]

    filteredTable(dt)
  })

  output$callTable <- DT::renderDataTable({

    req(filteredTable())
    return(DT::datatable(filteredTable(),
                         selection = "single",
                         style = "bootstrap",
                         class = DT:::DT2BSClass(c("compact", "hover", "stripe")),
                         options = list(columnDefs = list(list(
                          targets = seq_len(ncol(filteredTable())),
                          render = DT::JS(
                            "function(data, type, row, meta) {",
                            "if (data == null) return '';",
                            "return type === 'display' && data.length > 20 ?",
                            "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                            "}"
                          )))
                        )))
  })


  ### SIDEBAR

  output$numTotalRows    <- renderText({ paste("Total calls in table: ",    nrow(fullTable())) })
  output$numFilteredRows <- renderText({ paste("Filtered calls in table: ", nrow(filteredTable())) })

  ### REACTOME TAB

  observeEvent(input$openReactomeButton, {

    req(filteredTable())
    req(filteredTable()$symbol)

    hgncGenes <- unique(filteredTable()$symbol)
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

  observeEvent(input$callTable_rows_selected, {

    req(filteredTable())

    selectedRow <- filteredTable()[input$callTable_rows_selected]

    matchingGene <- shiny.huge.geneTable[
      symbol == selectedRow$Symbol |
      grepl(selectedRow$Symbol, prev_symbol) |
      grepl(selectedRow$Symbol, alias_symbol)
    ]

    showModal(modalDialog(
      shiny.huge.detailDivElement("Sample:", selectedRow$Sample),
      shiny.huge.detailDivElement("Chromosome:", selectedRow$Chr),
      shiny.huge.detailDivElement("Gene (found in table):", selectedRow$Symbol),
      shiny.huge.detailDivElement("Gene (current HGNC symbol):", matchingGene$symbol),
      shiny.huge.detailDivElement("Gene full name:", matchingGene$name),
      shiny.huge.detailDivElement("Gene description:", matchingGene$description),
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

  observeEvent(input$modalOkBtn, {
    removeModal()
  })
})
