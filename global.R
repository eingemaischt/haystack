loadOrInstall <- function (packageName, type = "CRAN") {

  isPackageInstalled <- packageName %in% rownames(installed.packages())

  if (!isPackageInstalled) {

    if (type == "CRAN") {

      install.packages(packageName)

    } else if (type == "bioc") {

      source("https://bioconductor.org/biocLite.R")
      biocLite(packageName)

    }

  }

  library(packageName, character.only = TRUE)
}

cranPackages <- c(
  "shiny",
  "shinydashboard",
  "data.table",
  "hashmap",
  "httr",
  "ggplot2",
  "openxlsx",
  "htmlwidgets",
  "DT",
  "VennDiagram",
  "xml2",
  "rclipboard"
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

source("annotation.R")