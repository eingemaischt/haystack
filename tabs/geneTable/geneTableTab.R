tabs.geneTable.module <- function (input, output, session, geneTable, filteredCallTable, expressionFilter) {

  callModule(tableDownload, "geneDownload", geneTable, "genes-")

  output$geneTable <- DT::renderDataTable({
    req(geneTable())

    return(DT::datatable(geneTable(),
                         selection = "single",
                         style = "bootstrap",
                         class = DT:::DT2BSClass(c("hover", "stripe")))
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
