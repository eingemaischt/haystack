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
        sliderInput("sampleNumber",  label = "Number of samples affected", min = 0, max = 100, value = c(0, 100), step = 1),
        sliderInput("minReadDepth", label = "Minimum read depth", min = 0, max = 200, value = 0),
        sliderInput("minVariantDepth", label = "Minimum variant depth", min = 0, max = 200, value = 0),
        sliderInput("readVariantFrequency", label = "Frequency (variant depth / read depth)", min = 0, max = 1, step = 0.05, value = c(0,1)),
        numericInput("maxAFPopmax", label = "Maximum AF Popmax", min = 0, max = 100, value = 100),
        checkboxGroupInput("genotypes", label = "Genotype", choices = list(
          "unknown" = "unknown",
          "0/0" = "hom_ref",
          "0/1, 1/0" = "het",
          "1/1" = "hom_alt"
        ), selected = c("unknown", "hom_ref", "het", "hom_alt")),
        checkboxInput("onlyCompoundHeterozygosity", label = "Show only compound heterozygosity candidates", value = FALSE),
        selectizeInput("expressions", label = "Tissue expression", multiple = TRUE, choices = unique(shiny.huge.gtexExpression$tissue)),
        sliderInput("scaledTPM", label = "Scaled Expression TPM value", min = 0, max = 1, value = c(0,1), step = 0.01),
        selectizeInput("consequences", label = "Consequences", multiple = TRUE, choices = "NA"),
        selectizeInput("studies", label = "Study", multiple = TRUE, choices = "NA"),
        actionButton("filterReset", label = "Reset filters"),
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
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
        tags$style(
          HTML('.wrapper {height: auto !important; position:relative; overflow-x:hidden; overflow-y:hidden}')
        ),
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
              width = 8,
              p("You can upload your variant call table here."),
              hr(),
              fileInput("callFile", "Choose CSV File",
                        accept = c(
                          "text/csv",
                          "text/comma-separated-values,text/plain",
                          ".csv")
              )
            ),
            box(
              title = "Download",
              width = 4,
              p("You can download the filtered variant call table here. Please note that only the sidebar filters apply to this file, not the sorting orders or search bar."),
              hr(),
              downloadButton("filteredCallTableDownload")
            )
          ),
          fluidRow(
            box(
              title = "Select displayed columns",
              solidHeader = TRUE,
              p("Select which columns to display here. Beware to not disable columns that you want to use otherwise (e.g. histograms)."),
              checkboxGroupInput("selectedColumns", label = "Show columns: ", inline = TRUE)
            ),
            box(
              width = 3,
              title = "Symbols with no expression",
              footer = "These symbols have a valid HGNC symbol but show no expression data in any tissue (or all tissues show tpm values of 0).",
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
              downloadButton("gtexDownload", label = "Download expression table (GTEx)")
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
              style = "overflow-y: scroll; max-height: 400px",
              title = "Symbols only in our data",
              verbatimTextOutput("ourSymbols"),
              uiOutput("ourSymbolsToClipboard")
            ),
            box(
              width = 4,
              style = "overflow-y: scroll; max-height: 400px",
              title = "Intersecting symbols",
              verbatimTextOutput("intersectingSymbols"),
              uiOutput("intersectingSymbolsToClipboard")
            ),
            box(
              width = 4,
              style = "overflow-y: scroll; max-height: 400px",
              title = "Symbols only in their data",
              verbatimTextOutput("theirSymbols"),
              uiOutput("theirSymbolsToClipboard")
            )
          ),
          fluidRow(
            box(
              width = 6,
              title = "Gene List Upload",
              footer = "You can upload a gene list here to compare to the genes found in the variant table. Every line must contain exactly one HGNC symbol.",
              fileInput("geneComparisonListUpload", label = "Upload genes for comparison")
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