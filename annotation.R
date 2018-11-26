library(data.table)
library(biomaRt)
library(rentrez)
library(hashmap)
library(hpar)

shiny.huge.martMirror <- "www.ensembl.org"
shiny.huge.allowedReliabilities <- c("Approved", "Supported", "Enhanced")

shiny.huge.geneTable <- (function () {

  hgncTable <- fread(cmd = "zcat annotation/hgnc_complete_set.txt")
  ncbiTable <- fread(cmd = "zcat annotation/ncbi-gene-descriptions.csv")

  geneTable <- merge(hgncTable, ncbiTable, by.x = "symbol", by.y = "hgnc_symbol", all = TRUE)

  approvedGeneTable <- geneTable[status == "Approved"]

  return(approvedGeneTable)
})()

shiny.huge.symbolToIndexMap <- (function () {

  symbolMap <- hashmap(shiny.huge.geneTable$symbol, seq_along(shiny.huge.geneTable$symbol))

  splitAliases <- strsplit(shiny.huge.geneTable$alias_symbol, "|", fixed = TRUE)
  splitPrevSymbols <- strsplit(shiny.huge.geneTable$prev_symbol, "|", fixed = TRUE)

  aliasIndices <- lapply(seq_along(splitAliases), function (i) rep(i, length(splitAliases[[i]])))
  prevIndices <- lapply(seq_along(splitPrevSymbols), function (i) rep(i, length(splitPrevSymbols[[i]])))

  splitAliases <- unlist(splitAliases)
  splitPrevSymbols <- unlist(splitPrevSymbols)
  aliasIndices <- unlist(aliasIndices)
  prevIndices <- unlist(prevIndices)

  secondaryNames <- c(splitAliases, splitPrevSymbols)
  secondaryIndices <- c(aliasIndices, prevIndices)

  validIndices <- !(duplicated(secondaryNames) | duplicated(secondaryNames, fromLast = TRUE)) & !secondaryNames %in% shiny.huge.geneTable$symbol

  secondaryKeys <- secondaryNames[validIndices]
  secondaryValues <- secondaryIndices[validIndices]

  symbolMap[[secondaryKeys]] <- secondaryValues

  containingLowercaseKeys <- symbolMap$keys()
  containingLowercaseKeys <- containingLowercaseKeys[toupper(containingLowercaseKeys) != containingLowercaseKeys]

  symbolMap[[toupper(containingLowercaseKeys)]] <- symbolMap[[containingLowercaseKeys]]

  return(symbolMap)
})()

shiny.huge.gtexExpression <- (function () {

  if (exists("shiny.huge.gtexExpression")) return(shiny.huge.gtexExpression)

  gtexFile <- "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct"

  rawGtexData <- fread(cmd = paste0("zcat annotation/", gtexFile))

  scaledGtexData <- rawGtexData
  scaledGtexData[, 3:ncol(scaledGtexData)] <- data.table(t(apply(scaledGtexData[,-(1:2)], 1, function (x) {

    maxValue <- max(x)

    if (maxValue == 0) return(rep(0, length(x)))

    return(x / maxValue)

  })))

  rawGtexData <- melt(rawGtexData, id.vars = c("gene_id", "Description"), variable.name = "tissue", value.name = "tpm")
  scaledGtexData <- melt(scaledGtexData, id.vars = c("gene_id", "Description"), variable.name = "tissue", value.name = "tpm_scaled")

  combined <- merge(rawGtexData, scaledGtexData)

  combined$symbol <- shiny.huge.geneTable$symbol[shiny.huge.symbolToIndexMap[[combined$Description]]]
  combined$tissue <- tolower(combined$tissue)

  return(combined[tpm > 0 & !is.na(symbol)])

})()


shiny.huge.hpaProteinExpession <- (function () {

  if (exists("shiny.huge.hpaProteinExpession")) return(shiny.huge.hpaProteinExpession)

  data("hpaNormalTissue")

  normalizedExpression <- data.table(hpaNormalTissue, stringsAsFactors = FALSE)

  hgncIndices <- shiny.huge.symbolToIndexMap[[normalizedExpression$Gene.name]]
  isValidIndex <- !is.na(hgncIndices)

  normalizedExpression <- normalizedExpression[isValidIndex]

  colnames(normalizedExpression) <- tolower(colnames(normalizedExpression))

  normalizedExpression$symbol <- shiny.huge.geneTable$symbol[hgncIndices[isValidIndex]]
  normalizedExpression$level <- factor(normalizedExpression$level, ordered = TRUE, levels = c("Not detected", "Low", "Medium", "High"))

  normalizedExpression <- normalizedExpression[reliability %in% shiny.huge.allowedReliabilities & level > "Not detected"]

  rm(hpaNormalTissue, envir = .GlobalEnv)
  return(normalizedExpression)

})()

shiny.huge.hpaRnaExpression <- (function () {

  if (exists("shiny.huge.hpaRnaExpression")) return(shiny.huge.hpaRnaExpression)

  data("rnaGeneTissue")

  normalizedExpression <- data.table(rnaGeneTissue, stringsAsFactors = FALSE)

  hgncIndices <- shiny.huge.symbolToIndexMap[[normalizedExpression$Gene.name]]
  isValidIndex <- !is.na(hgncIndices)

  normalizedExpression <- normalizedExpression[isValidIndex]

  normalizedExpression$tpm <- normalizedExpression$Value
  normalizedExpression$Value <- NULL

  normalizedExpression$tissue <- as.character(normalizedExpression$Sample)
  normalizedExpression$Sample <- NULL

  normalizedExpression$symbol <- shiny.huge.geneTable$symbol[hgncIndices[isValidIndex]]

  geneMaxTpm <- normalizedExpression[,list(maxTpm = max(tpm)),by = Gene.name]

  matchingIndices <- match(normalizedExpression$Gene.name, geneMaxTpm$Gene.name)

  scaledTpm <- normalizedExpression$tpm / geneMaxTpm$maxTpm[matchingIndices]

  normalizedExpression$tpm_scaled <- scaledTpm

  normalizedExpression <- normalizedExpression[tpm > 0]

  rm(rnaGeneTissue, envir = .GlobalEnv)
  return(normalizedExpression)

})()

shiny.huge.annotateFromEntrez <- function (hgncSymbols) {

  uniqueGenes <- unique(unlist(hgncSymbols))
  uniqueGenes <- uniqueGenes[!is.null(uniqueGenes)]

  ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl", host = shiny.huge.martMirror)
  biomartResponse <- getBM(attributes = c("hgnc_symbol", "entrezgene"), filters = "hgnc_symbol", values = uniqueGenes, mart = ensembl)

  biomartResponse <- biomartResponse[!is.na(biomartResponse$entrezgene) & !duplicated(biomartResponse$hgnc_symbol),]

  chunks <- split(biomartResponse$entrezgene, ceiling(seq_along(biomartResponse$entrezgene) / 100))

  entrezSummaries <- lapply(chunks, function (chunk) {

    count <- 0
    response <- NULL

    while (count < 100 && is.null(response)) {
      tryCatch({

        webHistory <- entrez_post(db = "gene", id = chunk)

        response <- entrez_summary(db = "gene", web_history = webHistory, version = "2.0", always_return_list = TRUE)
      }, error = function (e) {
        print(paste("error connecting to entrez, try nr. :", count))
      })

      count <- count + 1
    }

    if (count == 100) stop("Too many tries when getting rentrez data")

    return(response)
  })

  entrezSummaries <- as.character(unlist(lapply(entrezSummaries, function (data) lapply(data, function (esummary) esummary$summary))))

  summaries <- data.frame(
    gene = biomartResponse$hgnc_symbol,
    description = entrezSummaries,
    stringsAsFactors = FALSE
  )

  summaries <- summaries[nchar(summaries$gene) > 0,]

  return(merge(summaries, data.frame(gene = hgncSymbols)))
}