shinyServer(function(input, output, session) {

  callTableReactives <- callModule(tabs.callTable.module, "callTableTab")

  fullCallTable <- callTableReactives$fullCallTable
  filteredCallTable <- callTableReactives$filteredCallTable

  geneTable <- callModule(
    sidebar.module,
    "sidebarFiltering",
    fullCallTable,
    filteredCallTable,
    callTableReactives$selectedColumns,
    callTableReactives$expressionFilter
  )

  callModule(tabs.geneTable.module, "geneTableTab", geneTable, filteredCallTable, callTableReactives$expressionFilter)
  callModule(tabs.geneComparison.module, "geneComparisonTab", fullCallTable, geneTable)
  callModule(tabs.annotations.module, "annotationsTab")
  callModule(tabs.pathways.module, "pathwaysTab", geneTable)

})
