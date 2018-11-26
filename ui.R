shinyUI(

  dashboardPage(
    dashboardHeader(title = "Human Genetics Shiny Germ Cells"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Call table", tabName = "callTableTab", icon = icon("th-list")),
        menuItem("Genes", tabName = "geneTab", icon = icon("leaf")),
        menuItem("Gene comparison", tabName = "geneComparisonTab", icon = icon("refresh")),
        menuItem("Annotations", tabName = "annotationTab", icon = icon("globe")),
        menuItem("Pathways",   tabName = "pathwayTab",   icon = icon("random")),

        sidebarFilteringUI("sidebarFiltering")
      )
    ),
    dashboardBody(
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
        rclipboardSetup()
      ),
      tabItems(
        tabItem(tabName = "callTableTab",
          callTableTabUI("callTableTab")
        ),
        tabItem(tabName = "geneTab",
          geneTableTabUI("geneTableTab")
        ),
        tabItem(tabName = "annotationTab",
          annotationsTabUI("annotationsTab")
        ),
        tabItem(tabName = "geneComparisonTab",
          geneComparisonTabUI("geneComparisonTab")
        ),
        tabItem(tabName = "pathwayTab",
          pathwaysTabUI("pathwaysTab")
        )
      )
    )
  )
)