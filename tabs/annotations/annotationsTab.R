annotationsTab <- function (input, output, session) {

  callModule(tableDownload, "annotationDownload", function () shiny.huge.geneTable, "annotation-")
  callModule(tableDownload, "gtexDownload", function () shiny.huge.gtexExpression, "gtex-expression-")
  callModule(tableDownload, "hpaRnaDownload", function () shiny.huge.hpaRnaExpression, "hpa-rna-")
  callModule(tableDownload, "hpaProteinDownload", function () shiny.huge.hpaProteinExpession, "hpa-protein-")

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



}