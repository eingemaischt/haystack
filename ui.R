shinyUI(

  dashboardPage(
    dashboardHeader(title = "Human Genetics Shiny Germ Cells"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Call table", tabName = "callTableTab", icon = icon("th-list")),
        menuItem("Genes", tabName = "geneTab", icon = icon("barcode")),
        menuItem("Histograms", tabName = "histogramTab", icon = icon("signal ")),
        menuItem("Pathways",   tabName = "pathwayTab",   icon = icon("random")),
        numericInput("minSamplePercent",  label = "Maximum p-value", value = 0.05),
        selectizeInput("expressions", label = "Tissue expression", multiple = TRUE, choices = c("None")),
        hr(),
        fluidPage(
          fluidRow(
            column(width = 12, textOutput("numFilteredRows"))
          ),
          fluidRow(
            column(width = 12, textOutput("numTotalRows"))
          )
        )
      )
    ),
    dashboardBody(
      tabItems(
        tabItem(tabName = "callTableTab",
          fluidRow(
            box(
              style = "overflow-x: scroll",
              title = "Variant Call Table",
              width = 12,
              solidHeader = TRUE,
              collapsible = TRUE,
              DT::dataTableOutput("callTable")
            )
          ),
          fluidRow(
            box(
              title = "Upload",
              p("You can upload your variant call table here."),
              hr(),
              fileInput("callFile", "Choose CSV File",
                        accept = c(
                          "text/csv",
                          "text/comma-separated-values,text/plain",
                          ".csv")
              )
            )
          ),
          fluidRow(
            box(
              title = "Select displayed columns",
              solidHeader = TRUE,
              p("Select which columns to display here. Beware to not disable columns that you want to use otherwise (e.g. histograms)."),
              checkboxGroupInput("selectedColumns", label = "Show columns: ", inline = TRUE)
            )
          )
        ),
        tabItem(tabName = "pathwayTab",
          fluidRow(
            box(
              title = "Open in Reactome browser",
              solidHeader = TRUE,
              p("You can perform a pathway overrepresentation analysis here. Gene symbols are uploaded to Reactome's servers and the analysis website is opened in a new tab."),
              hr(),
              actionButton("openReactomeButton", label = "Open in  Reactome", class="btn-info btn-block"),
              uiOutput("reactomeOpeningScript")
            )
          )
        )
      )
    )
  )
)
