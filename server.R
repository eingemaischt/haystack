shinyServer(function(input, output, session) {

  fullCallTable <- reactiveVal()
  filteredCallTable <- reactiveVal()
  expressionFilter <- reactiveVal()
  selectedColumns <- reactiveVal()
  geneTable <- reactiveVal()

  callModule(tabs.callTable.module, "callTableTab", fullCallTable, filteredCallTable, expressionFilter, selectedColumns)
  callModule(sidebar.module, "sidebarFiltering", fullCallTable, filteredCallTable, selectedColumns, expressionFilter, geneTable)
  callModule(tabs.geneTable.module, "geneTableTab", geneTable, filteredCallTable, expressionFilter)
  callModule(tabs.geneComparison.module, "geneComparisonTab", fullCallTable, geneTable)
  callModule(tabs.annotations.module, "annotationsTab")
  callModule(tabs.pathways.module, "pathwaysTab", geneTable)

})
