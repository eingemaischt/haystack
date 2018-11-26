shinyServer(function(input, output, session) {

  callTableReactives <- callModule(callTable, "callTableTab")

  fullCallTable <- callTableReactives$fullCallTable
  filteredCallTable <- callTableReactives$filteredCallTable

  geneTable <- callModule(sidebarFiltering, "sidebarFiltering", fullCallTable, filteredCallTable, callTableReactives$selectedColumns)

  callModule(geneTableTab, "geneTableTab", geneTable, filteredCallTable)
  callModule(geneComparisonTab, "geneComparisonTab", fullCallTable, geneTable)
  callModule(annotationsTab, "annotationsTab")
  callModule(pathwaysTab, "pathwaysTab", geneTable)

})
