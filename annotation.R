library(data.table)
library(biomaRt)
library(rentrez)
library(hpar)

shiny.p7.martMirror <- "www.ensembl.org"
shiny.p7.allowedLevels <- c("High", "Medium", "Low")
shiny.p7.allowedReliabilities <- c("Approved", "Supported")
shiny.p7.minTPMLevel <- 0.01

data("hpaNormalTissue")
data("rnaGeneTissue")

shiny.p7.normalExpressions <- as.data.table(hpaNormalTissue[hpaNormalTissue$Level %in% shiny.p7.allowedLevels & hpaNormalTissue$Reliability %in% shiny.p7.allowedReliabilities,])
# order factors for later use
shiny.p7.normalExpressions$Level <- factor(shiny.p7.normalExpressions$Level, levels = c("Not detected", "Low", "Medium", "High"), ordered = TRUE)

shiny.p7.rnaExpressions <- as.data.table(rnaGeneTissue[rnaGeneTissue$Value >= shiny.p7.minTPMLevel,])

shiny.p7.gtexExpression <- (function() {

  if (exists("shiny.p7.gtexExpression")) return(shiny.p7.gtexExpression)

  gtexFile <- "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct"
  gtexData <- fread(gtexFile, stringsAsFactors = FALSE)
  gtexData <- melt(gtexData, id.vars = c("gene_id", "Description"), variable.name = "tissue", value.name = "tpm")
  gtexData$tissue <- tolower(gtexData$tissue)

  # convert "ENSG00000235373.1" to "ENSG00000235373"
  gtexData$gene_id <- sapply(strsplit(gtexData$gene_id, ".", fixed = TRUE), function (components) components[1])

  return(gtexData[tpm >= shiny.p7.minTPMLevel])
})()


shiny.p7.annotateFromDbSNP <- function (dt) {

  snp <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host = shiny.p7.martMirror)
  gene <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = shiny.p7.martMirror)

  snpAnnotation <- getBM( attributes=c("refsnp_id","minor_allele_freq", "set_name", "set_description", "allele"), filters="snp_filter", values=dt$RSID, mart = snp)

  ensemblIds <- getBM( attributes=c("refsnp_id", "ensembl_gene_stable_id"), filters="snp_filter", values=dt$RSID, mart = snp)
  ensemblIds <- ensemblIds[nchar(ensemblIds$ensembl_gene_stable_id) > 0,]

  hgncSymbols <- getBM( attributes = c( "ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = unique(ensemblIds$ensembl_gene_stable_id), mart = gene)
  hgncSymbols <- hgncSymbols[nchar(hgncSymbols$hgnc_symbol) > 0,]

  ensemblIds$hgnc_symbol <- hgncSymbols$hgnc_symbol[match(ensemblIds$ensembl_gene_stable_id, hgncSymbols$ensembl_gene_id, nomatch = NA)]

  mergedIds <- data.table(ensemblIds)[, list(
    ENSEMBL_GENES = paste(unique(ensembl_gene_stable_id[!is.na(ensembl_gene_stable_id)], collapse = ",")),
    HGNC_NAMES = paste(unique(hgnc_symbol[!is.na(hgnc_symbol)]), collapse = ",")
  ), by = refsnp_id]

  mergedAnnotations <- merge(snpAnnotation, mergedIds, all = TRUE, by = "refsnp_id")

  hitIndices <- match(mergedAnnotations$refsnp_id, dt$RSID)

  dt[hitIndices, "ALLELES"] <- mergedAnnotations$allele
  dt[hitIndices, "MAF"] <- mergedAnnotations$minor_allele_freq
  dt[hitIndices, "SET_NAME"] <- mergedAnnotations$set_name
  dt[hitIndices, "SET_DESCRIPTION"] <- mergedAnnotations$set_description

  dt[hitIndices, "HGNC_NAMES"] <- mergedAnnotations$HGNC_NAMES
  dt[hitIndices, "ENSEMBL_GENES"] <- mergedAnnotations$ENSEMBL_GENES

  # set hgnc and ensembl to "" instead of <NA> to simplify later filtering
  dt[is.na(HGNC_NAMES), HGNC_NAMES := ""]
  dt[is.na(ENSEMBL_GENES), ENSEMBL_GENES := ""]

  return(dt)
}

shiny.p7.annotateFromLocalDbFile <- function (dt, dbFile, idColumn, mafColumn, setName, setDescription) {

  localDb <- fread(dbFile, data.table = FALSE, stringsAsFactors = FALSE)

  foundVariants <- which(localDb[, idColumn] %in% dt$RSID)

  hitIndices <- match(localDb[foundVariants, idColumn], dt$RSID)

  dt[hitIndices, "MAF"] <- localDb[foundVariants, mafColumn]
  dt[hitIndices, "SET_NAME"] <- setName
  dt[hitIndices, "SET_DESCRIPTION"] <- setDescription

  return(dt)
}

shiny.p7.annotateFromEntrez <- function (dt) {

  splitGenes <- strsplit(dt$HGNC_NAMES, split = ",")

  uniqueGenes <- unique(unlist(splitGenes))
  uniqueGenes <- uniqueGenes[!is.null(uniqueGenes)]

  ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl", host = shiny.p7.martMirror)
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
    summary = entrezSummaries,
    stringsAsFactors = FALSE
  )

  summaries <- summaries[nchar(summaries$summary) > 0,]
  rownames(summaries) <- summaries$gene

  dt$GENE_DESCRIPTION <- sapply(splitGenes, function (splitGene) {

    uniqueGenes <- unique(splitGene[splitGene %in% summaries$gene])

    return(paste(uniqueGenes, summaries[uniqueGenes, "summary"], sep = ":", collapse = "   "))
  })

  return(dt)
}

shiny.p7.queryExpressionAnnotations <- function (ensemblStrings) {

  splitEnsemblIds <- strsplit(ensemblStrings, ",")

  expressions <- lapply(splitEnsemblIds, function (ids) {

    targetColnames <- c("tissue", "value")

    normalExpressions <- shiny.p7.normalExpressions[shiny.p7.normalExpressions$Gene %in% ids, c("Tissue", "Level")]
    rnaGeneExpressions <- shiny.p7.rnaExpressions[shiny.p7.rnaExpressions$Gene %in% ids, c("Sample", "Value")]
    gtexExpressions <- shiny.p7.gtexExpression[shiny.p7.gtexExpression$gene_id %in% ids, c("tissue", "tpm")]

    colnames(normalExpressions) <- targetColnames
    colnames(rnaGeneExpressions) <- targetColnames
    colnames(gtexExpressions) <- targetColnames

    return(list(
      normal = normalExpressions,
      rna    = rnaGeneExpressions,
      gtex   = gtexExpressions
    ))
  })

  return(expressions)
}

shiny.p7.annotateExpressions <- function (dt) {

  expressionsQueried <- shiny.p7.queryExpressionAnnotations(dt$ENSEMBL_GENES)

  joinedExpressions <- sapply(expressionsQueried, function (expressions) paste(unique(c(
    as.character(expressions$normal$tissue),
    as.character(expressions$rna$tissue),
    as.character(expressions$gtex$tissue)
  )), collapse = ","))

  dt$EXPRESSIONS <- joinedExpressions

  return(dt)

}

shiny.p7.replaceToRsIds <- function (dt) {

  mappingFile <- "IDMapping.txt"
  mappingTable <- fread(mappingFile, stringsAsFactors = FALSE, data.table = FALSE)
  mappingTable <- mappingTable[mappingTable$RsID != ".",]
  mappingTable$RsID <- sapply(strsplit(mappingTable$RsID, ","), function (ids) ids[1])
  rownames(mappingTable) <- mappingTable$Name

  matchingIndices <- dt$SNP %in% mappingTable$Name
  matchingIds <- dt$SNP[matchingIndices]

  rsIds <- dt$SNP
  rsIds[matchingIndices] <- mappingTable[matchingIds, "RsID"]

  dt$RSID <- rsIds

  return(dt)
}

shiny.p7.annotateEffect <- function (dt) {

  illuminaEffectAnnotationFile <- "InfiniumPsychArray-24v1-1_A1.annotated.txt"
  illuminaEffectAnnotation <- fread(illuminaEffectAnnotationFile, fill = TRUE)

  # this assumes all snps in dt are found in illumina data (which should be the case for our studies)
  matchIndices <- match(dt$SNP, illuminaEffectAnnotation$Name)

  dt$TRANSCRIPTS <- illuminaEffectAnnotation$`Transcript(s)`[matchIndices]
  dt$EFFECTS <- illuminaEffectAnnotation$`Mutation(s)`[matchIndices]

  return(dt)

}

shiny.p7.annotateSNPFile <- function (fileName, errorHandler) {

  dt <- fread(fileName)

  annotatedColumns <- c(
    "RSID",
    "ALLELES",
    "MAF",
    "TRANSCRIPTS",
    "EFFECTS",
    "SET_NAME",
    "SET_DESCRIPTION",
    "HGNC_NAMES",
    "GENE_DESCRIPTION",
    "ENSEMBL_GENES",
    "EXPRESSIONS",
    "IS_CANDIDATE"
  )

  if (all(annotatedColumns %in% colnames(dt))) return(dt)

  dt <- shiny.p7.replaceToRsIds(dt)

  dcmDBDescription <- "Dilated cardiomyopathy (DCM) is an important cause of heart failure with a strong familial component. We performed an exome-wide array-based association study (EWAS) to assess the contribution of missense variants to sporadic DCM: https://doi.org/10.1371/journal.pone.0172995"
  illuminaDbDescription <- "Proprietary data provided by Illumina"

  dt <- shiny.p7.annotateFromLocalDbFile(dt, "DCM-DB.txt", "SNP", "MAF", "DCM exome study", dcmDBDescription)
  dt <- shiny.p7.annotateFromLocalDbFile(dt, "ILLUMINA-DB.txt", "Name", "Minor Freq", "Illumina SNP data", illuminaDbDescription)

  dt <- shiny.p7.annotateFromDbSNP(dt)
  dt <- shiny.p7.annotateFromEntrez(dt)
  dt <- shiny.p7.annotateEffect(dt)
  dt <- shiny.p7.annotateExpressions(dt)

  dt$IS_CANDIDATE <- "none"

  return(dt)
}