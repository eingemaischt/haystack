loadOrInstall <- function (packageName, type = "CRAN") {

  isPackageInstalled <- packageName %in% rownames(installed.packages())

  if (!isPackageInstalled) {

    if (type == "CRAN") {

      install.packages(packageName)

    } else if (type == "bioc") {

      BiocManager::install(packageName)

    }

  }

  library(packageName, character.only = TRUE)
}

cranPackages <- c(
  "shiny",
  "shinydashboard",
  "shinyFeedback",
  "data.table",
  "hashmap",
  "httr",
  "ggplot2",
  "openxlsx",
  "htmlwidgets",
  "DT",
  "VennDiagram",
  "xml2",
  "rclipboard",
  "stringr",
  "BiocManager"
)

biocPackages <- c(
  "biomaRt",
  "rentrez",
  "hpar"
)

for (package in cranPackages) {

  loadOrInstall(package)

}

for (package in biocPackages) {

  loadOrInstall(package, type = "bioc")

}

# Prevent VennDiagram package from writing a log file each time it is called
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# make sure submissions always fit into the server
options(shiny.maxRequestSize=100*1024^2)

source("annotation/annotation.R")

source("util/errorModal.R")

source("util/copyToClipboardButton.R")

source("util/detailModal.R")

source("util/tableDownload.R")
source("util/tableDownloadUI.R")

source("tabs/callTable/callTableTab.R")
source("tabs/callTable/callTableTabUI.R")

source("tabs/geneTable/geneTableTab.R")
source("tabs/geneTable/geneTableTabUI.R")

source("tabs/geneComparison/geneComparisonTab.R")
source("tabs/geneComparison/geneComparisonTabUI.R")

source("tabs/annotations/annotationsTab.R")
source("tabs/annotations/annotationsTabUI.R")

source("tabs/pathways/pathwaysTab.R")
source("tabs/pathways/pathwaysTabUI.R")

source("tabs/help/helpTabUI.R")

source("sidebarFiltering/sidebarFiltering.R")
source("sidebarFiltering/sidebarFilteringUI.R")
