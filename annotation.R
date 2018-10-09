library(data.table)
library(biomaRt)
library(rentrez)
library(hashmap)
#library(hpar)

shiny.huge.martMirror <- "www.ensembl.org"
shiny.huge.allowedLevels <- c("High", "Medium", "Low")
shiny.huge.allowedReliabilities <- c("Approved", "Supported")

# data("hpaNormalTissue")
# data("rnaGeneTissue")

# shiny.huge.normalExpressions <- as.data.table(hpaNormalTissue[hpaNormalTissue$Level %in% shiny.huge.allowedLevels & hpaNormalTissue$Reliability %in% shiny.huge.allowedReliabilities,])
# order factors for later use
# shiny.huge.normalExpressions$Level <- factor(shiny.huge.normalExpressions$Level, levels = c("Not detected", "Low", "Medium", "High"), ordered = TRUE)

# shiny.huge.rnaExpressions <- as.data.table(rnaGeneTissue[rnaGeneTissue$Value >= shiny.huge.minTPMLevel,])

shiny.huge.loadGtexExpressions <- function (normalizationFn, globalVariableName, minTPM) {

  if (exists(globalVariableName)) return(get(globalVariableName))

  gtexFile <- "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct"
  gtexData <- fread(gtexFile, stringsAsFactors = FALSE)
  gtexData[, 3:ncol(gtexData)] <- data.table(t(apply(gtexData[,-(1:2)], 1, normalizationFn)))
  gtexData <- melt(gtexData, id.vars = c("gene_id", "Description"), variable.name = "tissue", value.name = "tpm")
  gtexData$tissue <- tolower(gtexData$tissue)

  # convert "ENSG00000235373.1" to "ENSG00000235373"
  gtexData$gene_id <- sapply(strsplit(gtexData$gene_id, ".", fixed = TRUE), function (components) components[1])

  return(gtexData)

}

shiny.huge.gtexExpressionRaw <- shiny.huge.loadGtexExpressions(function (x) x, "shiny.huge.gtexExpression")
shiny.huge.gtexExpressionScaled <- shiny.huge.loadGtexExpressions(scales::rescale, "shiny.huge.gtexExpressionScaled")

shiny.huge.geneTable <- (function () {

  hgncTable <- fread("hgnc_complete_set.txt")
  ncbiTable <- fread("ncbi-gene-descriptions.csv")

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

  return(symbolMap)
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

shiny.huge.queryExpressionAnnotations <- function (ensemblStrings) {

  expressions <- lapply(ensemblStrings, function (ids) {

    targetColnames <- c("tissue", "value")

    # normalExpressions <- shiny.huge.normalExpressions[shiny.huge.normalExpressions$Gene %in% ids, c("Tissue", "Level")]
    # rnaGeneExpressions <- shiny.huge.rnaExpressions[shiny.huge.rnaExpressions$Gene %in% ids, c("Sample", "Value")]
    gtexExpressions <- shiny.huge.gtexExpressionRaw[shiny.huge.gtexExpressionRaw$gene_id %in% ids, c("tissue", "tpm")]

    # colnames(normalExpressions) <- targetColnames
    # colnames(rnaGeneExpressions) <- targetColnames
    colnames(gtexExpressions) <- targetColnames

    return(list(
      # normal = normalExpressions,
      # rna    = rnaGeneExpressions,
      gtex   = gtexExpressions
    ))
  })

  return(expressions)
}