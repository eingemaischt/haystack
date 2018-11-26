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
              width = 6,
              p("You can upload your variant call table here. Make sure to upload only correctly formatted files (comma separated, periods as decimal seperators, variant table as first sheet for .xlsx files). Both .csv and .xlsx files are supported."),
              hr(),
              fileInput("callFile", "Choose file",
                        accept = c(
                          "text/csv",
                          "text/comma-separated-values,text/plain",
                          ".csv",
                          "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                          ".xlsx")
              )
            ),
            box(
              title = "Download",
              width = 6,
              style = "height: 100%",
              p("You can download the filtered variant call table here. Please note that only the sidebar filters apply to this file, not the sorting orders or search bar."),
              hr(),
              downloadButton("filteredCallTableCsvDownload", label = "Download as CSV"),
              downloadButton("filteredCallTableXlsxDownload", label = "Download as XLSX")
            )
          ),
          fluidRow(
            box(
              title = "Select displayed columns",
              p("Select which columns to display here. Beware to not disable columns that you want to use otherwise."),
              checkboxGroupInput("selectedColumns", label = "Show columns: ", inline = TRUE)
            ),
            box(
              width = 3,
              title = "Symbols with no expression",
              footer = "These symbols have a valid HGNC symbol but show no expression data in any tissue (or all tissues show tpm values of 0 / 'not detected' protein level).",
              tags$div(
                style = "overflow-y: scroll; max-height: 400px",
                verbatimTextOutput("unexpressedSymbols")
              ),
              tags$div(
                "Total: ",
                textOutput("totalUnexpressedSymbols", inline = TRUE)
              ),
              uiOutput("unexpressedToClipboard")
            ),
            box(
              width = 3,
              title = "Unrecognized symbols",
              footer = "These symbols are no valid (up-to-date) HGNC symbols. Also, no valid previous or alias name could be found.",
              tags$div(
                style = "overflow-y: scroll; max-height: 400px",
                verbatimTextOutput("unrecognizedSymbols")
              ),
              tags$div(
                "Total: ",
                textOutput("totalUnrecognizedSymbols", inline = TRUE)
              ),
              uiOutput("unrecognizedToClipboard")
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
              width = 4,
              solidHeader = TRUE,
              footer = "Please note that downloads always contain the whole tables. Filters and sorting orders do not apply.",
              title = "Download",
              downloadButton("annotationDownload", label = "Download gene annotation table (HGNC/NCBI)"),
              hr(),
              downloadButton("gtexDownload", label = "Download expression table (GTEx)"),
              hr(),
              downloadButton("hpaRnaDownload", label = "Download expression table (HPA RNA)"),
              hr(),
              downloadButton("hpaProteinDownload", label = "Download expression table (HPA protein)")
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