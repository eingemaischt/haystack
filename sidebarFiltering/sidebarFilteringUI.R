sidebarFilteringUI <- function (id) {

  ns <- NS(id)

  tagList(
    sliderInput(ns("sampleNumber"),  label = "Number of samples affected", min = 0, max = 100, value = c(0, 100), step = 1),
    sliderInput(ns("minReadDepth"), label = "Minimum read depth", min = 0, max = 200, value = 0),
    sliderInput(ns("minVariantDepth"), label = "Minimum variant depth", min = 0, max = 200, value = 0),
    sliderInput(ns("readVariantFrequency"), label = "Frequency (variant depth / read depth)", min = 0, max = 1, step = 0.05, value = c(0,1)),
    numericInput(ns("maxAFPopmax"), label = "Maximum AF Popmax", min = 0, max = 100, value = 100),
    checkboxGroupInput(ns("genotypes"), label = "Genotype", choices = list(
      "unknown" = "unknown",
      "0/1, 1/0" = "0/1, 1/0",
      "1/1" = "1/1"
    ), selected = c("unknown", "0/1, 1/0", "1/1")),
    checkboxInput(ns("onlyCompoundHeterozygosity"), label = "Show only compound heterozygosity candidates", value = FALSE),
    selectizeInput(ns("expressions"), label = "GTEx tissue expression", multiple = TRUE, choices = unique(
      shiny.huge.gtexExpression$tissue)),
    sliderInput(ns("scaledTPM"), label = "Scaled GTEx TPM value", min = 0, max = 1, value = c(0,1), step = 0.01),
    selectizeInput(ns("consequences"), label = "Consequence", multiple = TRUE, choices = "NA"),
    selectizeInput(ns("studies"), label = "Study", multiple = TRUE, choices = "NA"),
    selectizeInput(ns("chromosomes"), label = "Chromosome", multiple = TRUE, choices = "NA"),
    actionButton(ns("filterReset"), label = "Reset filters", icon = icon("repeat")),
    tableDownloadUI(ns("filterDownload"), "filters")
  )
}