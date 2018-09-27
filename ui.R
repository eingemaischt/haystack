shinyUI(

  dashboardPage(
    dashboardHeader(title = "Human Genetics Shiny Germ Cells"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Call table", tabName = "callTableTab", icon = icon("th-list")),
        menuItem("Genes", tabName = "geneTab", icon = icon("leaf")),
        menuItem("Annotations", tabName = "annotationTab", icon = icon("globe")),
        menuItem("Pathways",   tabName = "pathwayTab",   icon = icon("random")),
        sliderInput("sampleNumber",  label = "Number of samples affected", min = 0, max = 100, value = c(0, 100), step = 1),
        sliderInput("minReadDepth", label = "Minimum read depth", min = 0, max = 100, value = 0),
        sliderInput("minVariantDepth", label = "Minimum variant depth", min = 0, max = 100, value = 0),
        checkboxInput("onlyCompoundHeterozygosity", label = "Compound heterozygosity only", value = FALSE),
        selectizeInput("expressions", label = "Tissue expression", multiple = TRUE, choices = unique(shiny.huge.gtexExpression$tissue)),
        selectizeInput("consequences", label = "Consequences", multiple = TRUE, choices = "NA"),
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
        tabItem(tabName = "geneTab",
          fluidRow(
            box(
              style = "overflow-x: scroll",
              title = "Gene Table",
              width = 12,
              solidHeader = TRUE,
              collapsible = TRUE,
              DT::dataTableOutput("geneTable")
            )
          ),
          fluidRow(
            box(
              width = 6,
              solidHeader = TRUE,
              title = "Download",
              downloadButton("geneDownload", label = "Download table"),
              footer = "This download contains the filtered genes only."
            )
          )
        ),
        tabItem(tabName = "annotationTab",
          fluidRow(
            box(
              style = "overflow-x: scroll",
              title = "Annotation Table",
              width = 12,
              solidHeader = TRUE,
              footer = "This window shows the annotation data used for this app, it is extracted from the HGNC.",
              DT::dataTableOutput("annotationTable")
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
