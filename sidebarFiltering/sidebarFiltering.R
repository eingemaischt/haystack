sidebar.resetFilters <- function (session, callTable, geneFilter) {

  numberOfSamples <- length(unique(callTable$Sample))
  consequences <- unique(unlist(strsplit(callTable$Consequence, ",")))

  studies <- unique(callTable$Studie)

  minTpmRank <- min(annotation.gtexExpression$tpm_rank)
  maxTpmRank <- max(annotation.gtexExpression$tpm_rank)

  geneFilter("")

  maxGnomadOeScore <- max(callTable$`gnomAD oe`, na.rm = TRUE)
  maxGnomadOeScore <- ifelse(is.infinite(maxGnomadOeScore), 1, maxGnomadOeScore)

  updateSliderInput(session, "sampleNumber", value = c(0, numberOfSamples), max = numberOfSamples)
  updateSliderInput(session, "minReadDepth", value = 0)
  updateSliderInput(session, "minVariantDepth", value = 0)
  updateSliderInput(session, "readVariantFrequency", value = c(0,1))
  updateSliderInput(session, "balance", value = c(0, 100))
  updateSliderInput(session, "scaledTPM", value = c(0,1))
  updateSliderInput(session, "tpmRank", value = c(minTpmRank, maxTpmRank))
  updateSliderInput(session, "gnomadOe", value = c(0, maxGnomadOeScore), max = maxGnomadOeScore)
  updateTextAreaInput(session, "genes", value = "")
  updateSelectInput(session, "proteinLevel", selected = "Any")
  updateCheckboxGroupInput(session, "genotypes", selected = c("unknown", "0/1, 1/0", "1/1"))
  updateCheckboxInput(session, "onlyCompoundHeterozygosity", value = FALSE)
  updateNumericInput(session, "maxAFPopmax", value = 100)
  updateSelectizeInput(session, "expressions", selected = NULL, choices = unique(annotation.gtexExpression$tissue))
  updateSelectizeInput(session, "consequences", selected = NULL, choices = consequences)
  updateSelectizeInput(session, "studies", selected = NULL, choices = studies)
  updateSelectizeInput(session, "chromosomes", selected = NULL, choices = unique(callTable$Chr))
  updateSelectInput(session, "mpoPhenotypes", selected = NULL)

}

sidebar.intervalFilterString <- function(rangeInput) {

  return(paste0("[", rangeInput[1], ",", rangeInput[2], "]"))

}

sidebar.collapsedListFilter <- function (listFilter) {

  return(paste0(listFilter, collapse = ","))

}

sidebar.filterSettingsReactiveTable <- function (input, geneFilter) {

  return(function () {

    filters <- data.table(
      `Number of samples affected` = sidebar.intervalFilterString(input$sampleNumber),
      `Minimum read depth` = input$minReadDepth,
      `Minimum variant depth` = input$minVariantDepth,
      `Frequency (variant depth / read depth)` = sidebar.intervalFilterString(input$readVariantFrequency),
      `Balance` = sidebar.intervalFilterString(input$balance),
      `Maximum AF Popmax` = input$maxAFPopmax,
      `gnomAD observed/expected` = sidebar.intervalFilterString(input$gnomadOe),
      `Genotype` = sidebar.collapsedListFilter(input$genotypes),
      `Show only compound heterozygosity candidates` = input$onlyCompoundHeterozygosity,
      `GTEx tissue expression` = sidebar.collapsedListFilter(input$expressions),
      `MPO phenotypes` = sidebar.collapsedListFilter(input$mpoPhenotypes),
      `Scaled GTEx TPM value` = sidebar.intervalFilterString(input$scaledTPM),
      `GTEx TPM rank` = sidebar.intervalFilterString(input$tpmRank),
      `Gene` = sidebar.collapsedListFilter(geneFilter()),
      `Consequence` = sidebar.collapsedListFilter(input$consequences),
      `Study` = sidebar.collapsedListFilter(input$studies),
      `Chromosome` = sidebar.collapsedListFilter(input$chromosomes)
    )

    molten <- melt(filters, variable.name = "filter", measure.vars = colnames(filters))

    return(molten)
  })

}

sidebar.module <- function (input, output, session, fullCallTable, filteredCallTable, selectedColumns, expressionFilter, geneTable) {

  geneFilter <- reactiveVal("")

  callModule(util.tableDownload.module, "filterDownload", sidebar.filterSettingsReactiveTable(input, geneFilter), "filter-settings-")

  observeEvent(input$filterReset, {

    req(fullCallTable())

    sidebar.resetFilters(session, fullCallTable(), geneFilter)

  })

  observeEvent(input$maxAFPopmax, {

    # Input ID needs namespace prefix (unfortunately)
    feedbackDanger(
      inputId = "sidebarFiltering-maxAFPopmax",
      condition = is.na(input$maxAFPopmax),
      text = "Invalid number, make sure to use appropriate delimiter."
    )

  })

  observe({

    expressionFilter(input$expressions)

  })

  observe({

    req(fullCallTable())

    sidebar.resetFilters(session, fullCallTable(), geneFilter)

  })

  observe({
    req(fullCallTable())
    req(selectedColumns())
    req(input$sampleNumber)
    req(input$minReadDepth)
    req(input$minVariantDepth)
    req(input$maxAFPopmax)
    req(input$scaledTPM)
    req(input$tpmRank)
    req(input$gnomadOe)
    req(input$balance)

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
        (is.na(`gnomAD oe`) | is.nan(`gnomAD oe`) | (input$gnomadOe[1] <= `gnomAD oe` & input$gnomadOe[2] >= `gnomAD oe`)) &
        (is.na(`Balance`) | (input$balance[1] <= `Balance` & input$balance[2] >= `Balance`)) &
        (is.null(input$consequences) | grepl(paste0(input$consequences, collapse = "|"), Consequence)) &
        (is.null(input$studies) | grepl(paste0(input$studies, collapse = "|"), Studie)) &
        (is.null(input$chromosomes) | Chr %in% input$chromosomes),
      selectedColumns(), with = FALSE
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

    matchingGeneIndices <- annotation.symbolToIndexMap[[gt$Symbol]]
    matchingGeneNames <- annotation.geneTable$symbol[matchingGeneIndices]

    gt$name <- annotation.geneTable$name[matchingGeneIndices]
    gt$locus_type <- annotation.geneTable$locus_type[matchingGeneIndices]
    gt$family <- annotation.geneTable$gene_family[matchingGeneIndices]
    gt$description <- annotation.geneTable$description[matchingGeneIndices]

    gtexFilteredGenes <- annotation.gtexExpression[
      tpm_scaled >= input$scaledTPM[1] & tpm_scaled <= input$scaledTPM[2] &
      tpm_rank >= input$tpmRank[1] & tpm_rank <= input$tpmRank[2] &
        tissue %in% input$expressions
    ]

    genesKept <- unlist(strsplit(input$genes, "\n"))
    genesKept <- str_replace_all(genesKept, " ", "")
    genesKept <- genesKept[nchar(genesKept) > 0]

    geneFilter(genesKept)

    genesKeptIndices <- annotation.symbolToIndexMap[[genesKept]]
    genesKeptValidNames <- !is.na(genesKeptIndices)
    genesKeptNames <- annotation.geneTable$symbol[genesKeptIndices[genesKeptValidNames]]

    genesKeptInvalidNames <- toupper(genesKept[!genesKeptValidNames])

    mpoPhenotypes <- annotation.mpoPhenotypes[symbol %in% matchingGeneNames]

    groupedPhenotypes <- mpoPhenotypes[, list(phenotype = paste0(target, collapse = ", ")), by = symbol]
    mpoIndices <- match(matchingGeneNames, groupedPhenotypes$symbol)

    gt$mpoPhenotypes <- groupedPhenotypes$phenotype[mpoIndices]

    mpoFilteredGenes <- mpoPhenotypes[target %in% input$mpoPhenotypes]$symbol

    gt <- gt[
      input$sampleNumber[1] <= samples &
      input$sampleNumber[2] >= samples &
      (is.null(input$expressions) | matchingGeneNames %in% gtexFilteredGenes$symbol) &
      (is.null(input$mpoPhenotypes) | matchingGeneNames %in% mpoFilteredGenes) &
      (length(genesKept) == 0 | (!is.na(matchingGeneIndices) & matchingGeneNames %in% genesKeptNames) | (is.na(matchingGeneIndices) & toupper(Symbol) %in% genesKeptInvalidNames))
    ]
    ct <- ct[Symbol %in% gt$Symbol]

    geneTable(gt)
    filteredCallTable(ct)

  })


}
