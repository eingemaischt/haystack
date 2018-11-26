resetFilters <- function (session, callTable) {

  numberOfSamples <- length(unique(callTable$Sample))
  consequences <- unique(unlist(strsplit(callTable$Consequence, ",")))

  studies <- unique(callTable$Studie)

  # TODO: column selection?
  #updateCheckboxGroupInput(session, "selectedColumns", choices = colnames(callTable), selected = colnames(callTable))
  updateSliderInput(session, "sampleNumber", value = c(0, numberOfSamples), max = numberOfSamples)
  updateSliderInput(session, "minReadDepth", value = 0)
  updateSliderInput(session, "minVariantDepth", value = 0)
  updateSliderInput(session, "readVariantFrequency", value = c(0,1))
  updateSliderInput(session, "scaledTPM", value = c(0,1))
  updateSelectInput(session, "proteinLevel", selected = "Any")
  updateCheckboxGroupInput(session, "genotypes", selected = c("unknown", "0/1, 1/0", "1/1"))
  updateCheckboxInput(session, "onlyCompoundHeterozygosity", value = FALSE)
  updateNumericInput(session, "maxAFPopmax", value = 100)
  updateSelectizeInput(session, "expressions", selected = NULL, choices = unique(shiny.huge.gtexExpression$tissue))
  updateSelectizeInput(session, "consequences", selected = NULL, choices = consequences)
  updateSelectizeInput(session, "studies", selected = NULL, choices = studies)
  updateSelectizeInput(session, "chromosomes", selected = NULL, choices = unique(callTable$Chr))

}

intervalFilterString <- function(rangeInput) {

  return(paste0("[", rangeInput[1], ",", rangeInput[2], "]"))

}

collapsedListFilter <- function (listFilter) {

  return(paste0(listFilter, collapse = ","))

}

filterSettingsReactiveTable <- function (input) {

  return(function () {

    filters <- data.table(
      `Number of samples affected` = intervalFilterString(input$sampleNumber),
      `Minimum read depth` = input$minReadDepth,
      `Minimum variant depth` = input$minVariantDepth,
      `Frequency (variant depth / read depth)` = intervalFilterString(input$readVariantFrequency),
      `Maximum AF Popmax` = input$maxAFPopmax,
      `Genotype` = collapsedListFilter(input$genotypes),
      `Show only compound heterozygosity candidates` = input$onlyCompoundHeterozygosity,
      `GTEx tissue expression` = collapsedListFilter(input$expressions),
      `Scaled GTEx TPM value` = intervalFilterString(input$scaledTPM),
      `Consequence` = collapsedListFilter(input$consequences),
      `Study` = collapsedListFilter(input$studies),
      `Chromosome` = collapsedListFilter(input$chromosomes)
    )

    molten <- melt(filters, variable.name = "filter", measure.vars = colnames(filters))

    return(molten)
  })

}

sidebarFiltering <- function (input, output, session, fullCallTable, selectedColumns) {


  filteredCallTable <- reactiveVal()
  geneTable <- reactiveVal()

  observeEvent(input$filterReset, {

    req(fullCallTable())

    resetFilters(session, fullCallTable())

  })

  # TODO: fix table download
  #output$filterDownload <- shiny.huge.handleTableDownload(shiny.huge.filterSettingsReactiveTable(input), "filter-settings-")

  observe({

    maxPopMax <- input$maxAFPopmax

    # TODO: handle error notification
    #shiny.huge.handleErrorNotification(maxPopMax, "AF Popmax ", "popmaxErrorNotification")

  })

  observe({

    req(fullCallTable())

    resetFilters(session, fullCallTable())

  })

  observe({

    req(fullCallTable())
    # TODO: column selection?
    #req(selectedColumns)
    req(input$sampleNumber)
    req(input$minReadDepth)
    req(input$minVariantDepth)
    req(input$maxAFPopmax)
    req(input$scaledTPM)

    ct <- fullCallTable()

    variantFrequency <- ct$`Variant depth` / ct$`Read depth`

    ct <- ct[
      (
        ("0/1, 1/0" %in% input$genotypes & (Genotype == "0/1" | Genotype == "1/0")) |
          ("1/1" %in% input$genotypes & Genotype == "1/1") |
          ("unknown" %in% input$genotypes & (grepl(".", Genotype, fixed = TRUE) | Genotype == "0/0"))
      ) &
        `Read depth` >= input$minReadDepth &
        (`Variant depth` >= input$minVariantDepth | is.na(`Variant depth`)) &
        ((input$readVariantFrequency[1] <= variantFrequency & variantFrequency <= input$readVariantFrequency[2]) | `Read depth` == 0) &
        (is.na(`AF Popmax`) | input$maxAFPopmax >= `AF Popmax`) &
        (is.null(input$consequences) | grepl(paste0(input$consequences, collapse = "|"), Consequence)) &
        (is.null(input$studies) | grepl(paste0(input$studies, collapse = "|"), Studie)) &
        (is.null(input$chromosomes) | Chr %in% input$chromosomes)#,
      # TODO: column selection?
      #selectedColumns, with = FALSE
      ]

    # Filter for compound heterozygosity only after other variant filters have been applied
    # Otherwise, there will be variants that are not compound based on the filters
    sampleSymbolStrings <- paste0(ct$Sample, ct$Symbol)
    sampleSymbolDuplicates <- duplicated(sampleSymbolStrings) | duplicated(sampleSymbolStrings, fromLast = TRUE)

    ct <- ct[!input$onlyCompoundHeterozygosity | sampleSymbolDuplicates]

    # count mutations per gene per sample by multiple aggregations
    gt <- ct[, list(samples = .N), by = .(Symbol, Sample)]
    gt <- gt[, list(samples = .N, compound_het_samples = sum(samples > 1)), by = .(Symbol)]
    gt <- gt[order(samples, decreasing = TRUE)]

    matchingGeneIndices <- shiny.huge.symbolToIndexMap[[gt$Symbol]]
    matchingGeneNames <- shiny.huge.geneTable$symbol[matchingGeneIndices]

    gt$name <- shiny.huge.geneTable$name[matchingGeneIndices]
    gt$locus_type <- shiny.huge.geneTable$locus_type[matchingGeneIndices]
    gt$family <- shiny.huge.geneTable$gene_family[matchingGeneIndices]
    gt$description <- shiny.huge.geneTable$description[matchingGeneIndices]

    gtexFilteredGenes <- shiny.huge.gtexExpression[
      tpm_scaled >= input$scaledTPM[1] & tpm_scaled <= input$scaledTPM[2] &
        tissue %in% input$expressions
      ]

    gt <- gt[
      input$sampleNumber[1] <= samples &
        input$sampleNumber[2] >= samples &
        (is.null(input$expressions) | matchingGeneNames %in% gtexFilteredGenes$symbol)
      ]
    ct <- ct[Symbol %in% gt$Symbol]

    geneTable(gt)
    filteredCallTable(ct)

  })

  return(list(
    geneTable = geneTable,
    filteredCallTable = filteredCallTable
  ))
}