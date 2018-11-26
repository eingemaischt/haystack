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
              tableDownloadUI("geneDownload", "genes"),
              footer = "This download contains the filtered genes only. Please note that only the sidebar filters apply to this file, not the sorting orders or search bar."
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
          ),
          fluidRow(
            box(
              width = 6,
              solidHeader = TRUE,
              footer = "Please note that downloads always contain the whole tables. Filters and sorting orders do not apply.",
              title = "Download",
              tableDownloadUI("annotationDownload", "gene annotations"),
              hr(),
              tableDownloadUI("gtexDownload", "GTEx table"),
              hr(),
              tableDownloadUI("hpaRnaDownload", "HPA RNA table"),
              hr(),
              tableDownloadUI("hpaProteinDownload", "HPA protein table")
            )
          )
        ),
        tabItem(tabName = "geneComparisonTab",
          fluidRow(
            box(
              width = 12,
              title = "Venn Diagram",
              footer = "You can save this diagram using right click -> 'Save image as'.",
              plotOutput("geneComparisonVenn")
            )
          ),
          fluidRow(
            box(
              width = 4,
              title = "Symbols only in our data",
              tags$div(
                style = "overflow-y: scroll; max-height: 400px",
                verbatimTextOutput("ourSymbols")
              ),
              hr(),
              uiOutput("ourSymbolsToClipboard")
            ),
            box(
              width = 4,
              title = "Intersecting symbols",
              tags$div(
                style = "overflow-y: scroll; max-height: 400px",
                verbatimTextOutput("intersectingSymbols")
              ),
              hr(),
              uiOutput("intersectingSymbolsToClipboard")
            ),
            box(
              width = 4,
              title = "Symbols only in their data",
              tags$div(
                style = "overflow-y: scroll; max-height: 400px",
                verbatimTextOutput("theirSymbols")
              ),
              hr(),
              uiOutput("theirSymbolsToClipboard")
            )
          ),
          fluidRow(
            box(
              width = 4,
              title = "Gene List Upload",
              footer = "You can upload a gene list here to compare to the genes found in the variant table. Every line must contain exactly one HGNC symbol.",
              fileInput("geneComparisonListUpload", label = "Upload genes for comparison")
            ),
            box(
              width = 3,
              title = "Filters",
              footer = "Choose here if sidebar filters should be applied to symbols in our data or if all symbols should be considered regardless of filtering.",
              checkboxInput("geneComparisonUseFilters", value = TRUE, label = "Apply sidebar filters to our symbols")
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