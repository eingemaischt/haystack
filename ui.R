shinyUI(

  dashboardPage(
    dashboardHeader(title = "Haystack"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Call table", tabName = "callTableTab", icon = icon("th-list")),
        menuItem("Genes", tabName = "geneTab", icon = icon("leaf")),
        menuItem("Gene comparison", tabName = "geneComparisonTab", icon = icon("refresh")),
        menuItem("Annotations", tabName = "annotationTab", icon = icon("globe")),
        menuItem("Pathways",   tabName = "pathwayTab",   icon = icon("random")),
        menuItem("Help", tabName = "helpTab", icon = icon("question-circle")),

        sidebar.moduleUI("sidebarFiltering")
      )
    ),
    dashboardBody(
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
        useShinyFeedback(),
        rclipboardSetup()
      ),
      tabItems(
        tabItem(tabName = "callTableTab",
          tabs.callTable.moduleUI("callTableTab")
        ),
        tabItem(tabName = "geneTab",
          tabs.geneTable.moduleUI("geneTableTab")
        ),
        tabItem(tabName = "annotationTab",
          tabs.annotations.moduleUI("annotationsTab")
        ),
        tabItem(tabName = "geneComparisonTab",
          tabs.geneComparison.moduleUI("geneComparisonTab")
        ),
        tabItem(tabName = "pathwayTab",
          tabs.pathways.moduleUI("pathwaysTab")
        ),
        tabItem(tabName = "helpTab",
          tabs.help.moduleUI("helpTab")
        )
      )
    )
  )
)
