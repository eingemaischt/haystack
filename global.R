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
  "htmlwidgets",
  "DT",
  "xml2"
)

biocPackages <- c(
  "biomaRt",
  "rentrez"
  #"hpar"
)

for (package in cranPackages) {

  loadOrInstall(package)

}

for (package in biocPackages) {

  loadOrInstall(package, type = "bioc")

}

source("annotation.R")