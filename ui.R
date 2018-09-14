shinyUI(

  dashboardPage(
    dashboardHeader(title = "Male Germ Cells P7"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("SNP table",  tabName = "snpTableTab",   icon = icon("th-list")),
        menuItem("Histograms", tabName = "histogramTab", icon = icon("signal ")),
        menuItem("Pathways",   tabName = "pathwayTab",   icon = icon("random")),
        numericInput("maxPValue",  label = "Maximum p-value", value = 0.05),
        numericInput("minAbsBeta", label = "Minimum absolute beta", value = 0),
        numericInput("minNMiss", label = "Minimum non-missing genotypes", value = 0),
        numericInput("minR2", label = "Minimum R2", value = 0),
        selectizeInput("expressions", label = "Tissue expression", multiple = TRUE, choices = c("None")),
        selectInput("candidateFilter", label = "Show candidates", multiple = TRUE, choices = c("approved", "discarded", "none")),
        radioButtons("emptyGenesRadio", choices = c(include = TRUE, exclude = FALSE), label = "SNPs with no genes"),
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
        tabItem(tabName = "snpTableTab",
          fluidRow(
            box(
              style = "overflow-x: scroll",
              title = "SNP data table",
              width = 12,
              solidHeader = TRUE,
              collapsible = TRUE,
              DT::dataTableOutput("snpTable")
            )
          ),
          fluidRow(
            box(
              title = "Upload",
              p("You can upload your SNP table here. If the SNP table already contains all annotation columns, database queries are skipped."),
              hr(),
              fileInput("snpFile", "Choose CSV File",
                        accept = c(
                          "text/csv",
                          "text/comma-separated-values,text/plain",
                          ".csv")
              )
            ),
            box(
              title = "Download",
              p("You can download the filtered SNP table here. The downloaded file also contains database annotations. Whenever you re-upload a SNP table that already contains all annotated columns, no database needs to be queried, thus speeding up the loading process."),
              hr(),
              downloadButton("downloadAllData", label = "Download whole SNP table"),
              downloadButton("downloadFilteredData", label = "Download filtered SNP table")
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
        tabItem(tabName = "histogramTab",
          fluidRow(
            box(
              title = "BETA",
              solidHeader = TRUE,
              plotOutput("betaHistPlot"),
              p("BETA denotes the slope of the linear regression for the continuous output variable (e.g. FSH levels) dependent of the sample's genotype (e.g. GG, GT and TT). The bigger the BETA value, the higher the difference of the output variable (e.g. FSH) for different genotypes.")
            ),
            box(
              title = "p-value",
              solidHeader = TRUE,
              plotOutput("pHistPlot"),
              p("P is the p-value based on a Wald test during the quantitive association analysis. The smaller this value, the smaller the chance of randomly observing the relation between output variable and genotype.")
            )
          ),
          fluidRow(
            box(
              title = "R2",
              solidHeader = TRUE,
              plotOutput("r2Plot"),
              p("R2 (r squared) is the coefficient of determination that describes the variance of the output variable (e.g. FSH levels) explained by the genotype. Its value is between 0 and 1, where 1 means all variance may be explained by the SNP. The bigger r², the more variance explained by the SNP.")
            ),
            box(
              title = "NMISS",
              solidHeader = TRUE,
              plotOutput("nmissPlot"),
              p("NMISS describes the number of non-missing genotypes at each SNP. Genotypes might for example be missing because of too low signals at a SNP probe.")
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
