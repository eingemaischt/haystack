tabs.geneTable.module <- function (input, output, session, geneTable, filteredCallTable, expressionFilter) {

  callModule(util.tableDownload.module, "geneDownload", geneTable, "genes-")

  observe({

    req(geneTable())

    geneNames <- geneTable()$Symbol

    output$geneNamesToClipboard <- util.copyToClipboardButton(paste0(geneNames, collapse = "\n"), "geneNamesToClipboardBtn")

  })

  output$geneTable <- DT::renderDataTable({
    req(geneTable())

    return(DT::datatable(geneTable(),
                         selection = "single",
                         style = "bootstrap",
                         class = DT:::DT2BSClass(c("hover", "stripe")),
                         options = list(
                           columnDefs =
                             list(
                               list(
                                 visible = FALSE,
                                 targets = 8
                               )
                             )
                         ))
    )
  })

  observeEvent(input$geneTable_rows_selected,
               util.detailModal.showDetailModal(
                 geneTable()[input$geneTable_rows_selected, Symbol],
                 filteredCallTable,
                 expressionFilter(),
                 output,
                 "geneTableTab")
  )

  observeEvent(input$modalOkBtn, {
    removeModal()
  })

}
