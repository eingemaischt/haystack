tabs.annotations.module <- function (input, output, session) {

  callModule(tableDownload, "annotationDownload", function () annotation.geneTable, "annotation-")
  callModule(tableDownload, "gtexDownload", function () annotation.gtexExpression, "gtex-expression-")
  callModule(tableDownload, "hpaRnaDownload", function () annotation.hpaRnaExpression, "hpa-rna-")
  callModule(tableDownload, "hpaProteinDownload", function () annotation.hpaProteinExpession, "hpa-protein-")

  ### ANNOTATION TABLE TAB

  output$annotationTable <- DT::renderDataTable(DT::datatable(annotation.geneTable,
                                                              selection = "single",
                                                              style = "bootstrap",
                                                              class = DT:::DT2BSClass(c("compact", "hover", "stripe")),
                                                              options = list(columnDefs = list(list(
                                                                targets = seq_len(ncol(annotation.geneTable)),
                                                                render = DT::JS(
                                                                  "function(data, type, row, meta) {",
                                                                  "if (data == null) return '';",
                                                                  "return type === 'display' && data.length > 20 ?",
                                                                  "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                                                                  "}"
                                                                )))
                                                              )))



}
