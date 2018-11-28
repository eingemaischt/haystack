tabs.help.moduleUI <- function (id) {

  ns <- NS(id)

  tagList(
    fluidRow(
      box(
        title = "Uploading Files",
        width = 4,
        tags$div("You can upload your Sciobase variant tables using the 'Call Table' tab. Both CSV and XLSX files are allowed. All uploaded tables are sanity-checked prior to loading. This might lead to potential errors when uploading."),
        tags$div("If you are uploading a file, make sure your table meets the following criteria:"),
        tags$ul(
          tags$li("Sciobase category information (usually the second row) must be removed prior to uploading."),
          tags$li("Sorting order should be patient-wise."),
          tags$li("Feature ID should be picked by variant effect predictor."),
          tags$li("There should be no samples where the sampe patient has been sequenced multiple times (e.g. both WES and WGS)."),
          tags$li("Dots must be used as decimal seperators. This is not done by default in some Excel configuration and may be changed using Excel options. This only applies for CSV files."),
          tags$li("For XLSX files with multiple sheets, make sure the variant table is the first sheet in the file so the server can recognize it.")
        )
      ),
      box(
        title = "Tabs",
        width = 4,
        tags$div("This website provides six tabs that you can use for gene hunting (located on the upper left hand side):"),
        tags$ul(
          tags$li(
            "Call table",
            tags$ul(
              tags$li("Allows upload of Sciobase variant tables. Variants may be sorted by clicking the respective column headers. Also, a search bar is provided on the upper right corner. Please note that this search bar has no effect on other parts of the website."),
              tags$li("Unrecognized/Unexpressed gene symbols are shown in the lower right corner."),
              tags$li("Variant call columns may optionally be hidden. Please make sure to hide no corners that are required for internal analyses (such as the Symbol column)."),
              tags$li("You can also download your filtered variant table again for further use. Be careful to import CSV files through Excel's user interface, as there may be conversion errors otherwise."),
              tags$li("You can click on every row in the table to get a detail view of the gene of the selected variant.")
            )
          ),
          tags$li(
            "Genes",
            tags$ul(
              tags$li("Gives a comprehensive overview ofer the genes from the variant table."),
              tags$li("Contains gene descriptions and locations for each gene from the variant table."),
              tags$li("Like with variant table rows, you may click on each gene table row to get a detailed overview over the gene selected.")
            )
          ),
          tags$li(
            "Gene comparison",
            tags$ul(
              tags$li("Here you can compare the genes from your data with gene lists from other sources."),
              tags$li("Make sure the second file contains only a list of gene names, and NOT a full variant table.")
            )
          ),
          tags$li(
            "Annotations",
            tags$ul(
              tags$li("Here you can download all tables that are used by the server for annotation and visualization purposes.")
            )
          ),
          tags$li(
            "Pathways",
            tags$ul(
              tags$li("Compare the genes from the filtered gene list to pathways using a ReactomeDB overrepresentationanalysis."),
              tags$li("Redirects to the ReactomeDB website, so make sure you allow popups when using this feature.")
            )
          ),
          tags$li(
            "Help",
            tags$ul(
              tags$li("Contains the information you are currently reading =)")
            )
          )
        )
      ),
      box(
        title = "Filtering",
        width = 4,
        tags$div("The main feature of this website consists of filtering the variant tables based on several criteria. Please note that the following rules apply for filtering:"),
        tags$ul(
          tags$li("All sidebar filters affect all tabs. So whenever you are setting filters for minimum read depth, the gene table will for example contain only those genes that remain after filtering."),
          tags$li("There are two ways to filter for expression data: using scaled TPM values and using ranked TPM values. The scaled values are computed by dividing TPM values for all tissues by the maximum TPM value in one tissue, so all values lie between 0 and 1. The ranked TPM values are simply computed by determining the order of the raw TPM values."),
          tags$li("TPM filters only apply if a tissue is selected using 'GTEx tissue expression'. If multiple tissues are selected, valid values in any of these tissues are sufficient for a gene to be kept."),
          tags$li("Expression filtering is only based on GTEx data, even if the detail plot also shows HPA data."),
          tags$li("All sample-based filtering is done using the Sample column, not the column for PatientNr's."),
          tags$li("Unknown genotype filters include genotypes such as 0/0, 0/., ./1 etc."),
          tags$li("Filtering by compound-heterozygosity often requires manual inspection: Sometimes, two variants lie exactly next to each other, so an artefact may be likely."),
          tags$li("Genotypes are extracted from the input variant table without checking for sex chromosomes, so heterozygosity and homozygosity might not work as expected for these variants.")
        )
      )
    )
  )

}