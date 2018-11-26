detailDivElement <- function (label, value) {

  return(tags$div(
    tags$b(label),
    tags$em(value)
  ))

}

modalExpressionPlot <- function (expressions, expressionFilter, title) {

  ### WORKAROUND
  ### when no expressions are found, just render some dummy plot so no errors are thrown
  ### see https://stackoverflow.com/questions/19918985/r-plot-only-text
  if (nrow(expressions) == 0) {

    errorMessage <- "No information available"

    return(renderPlot(
      ggplot() + annotate("text", x = 4, y = 25, size = 8, label = errorMessage)+
        theme_bw() +
        theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
    ))
  }
  ### END WORKAROUND

  return(renderPlot(
    ggplot(expressions, aes(x = reorder(tissue, value, FUN = max), y = value, fill = tissue %in% expressionFilter, colour = tissue %in% expressionFilter)) +
      geom_col() +
      coord_flip() +
      ggtitle(title) +
      xlab("Tissue") +
      ylab("Expression value") +
      scale_fill_manual(values = c("TRUE" = "springgreen", "FALSE" = "steelblue")) +
      scale_colour_manual(values = c("TRUE" = "springgreen4", "FALSE" = "steelblue4")) +
      guides(fill = FALSE, colour = FALSE),
    res = 90
  )
  )

}

geneExpressionModal <- function (selectedSymbol, callTableReactiveVal, expressionFilter, output, id) {

  ns <- NS(id)

  return({

    req(callTableReactiveVal())

    matchingGene <- shiny.huge.geneTable[shiny.huge.symbolToIndexMap[[selectedSymbol]]]

    showModal(modalDialog(
      detailDivElement("Location:", matchingGene$location),
      detailDivElement("Gene (found in table):", selectedSymbol),
      detailDivElement("Gene (current HGNC symbol):", matchingGene$symbol),
      detailDivElement("Gene full name:", matchingGene$name),
      detailDivElement("Gene family:", matchingGene$gene_family),
      detailDivElement("Gene description:", matchingGene$description),
      detailDivElement("Location type:", matchingGene$locus_type),
      detailDivElement("Ensembl ID:", matchingGene$ensembl_gene_id),
      tags$hr(),
      tags$b("Expression:"),
      tabsetPanel(
        tabPanel("GTex (TPM scaled)", plotOutput(ns("modalGTExScaledExpression"), height = "640px")),
        tabPanel("HPA RNA (TPM scaled)", plotOutput(ns("modalHpaRnaScaledExpression"), height = "640px")),
        tabPanel("HPA Protein", plotOutput(ns("modalHpaProteinExpression"), height = "640px")),
        tabPanel("Mutation types", tags$div(
          style = "overflow-x: scroll",
          tags$h5(paste0("Mutations for ", selectedSymbol), ":"),
          tableOutput(ns("modalMutationTypes"))
        )
        )
      ),
      title = "Details",
      footer = actionButton(ns("modalOkBtn"), label = "OK", icon = icon("ok")),
      size = "l",
      easyClose = TRUE
    ))

    callsInGene <- callTableReactiveVal()[Symbol == selectedSymbol]

    uniqueMutations <- callsInGene[,
                                   list("Samples" = length(unique(Sample))),
                                   by = list(HGVSc, Chr, Position, Consequence, `AF Popmax`)
                                   ]

    matchingGtexExpression <- shiny.huge.gtexExpression[symbol == matchingGene$symbol]
    matchingHpaRnaExpression <- shiny.huge.hpaRnaExpression[symbol == matchingGene$symbol]
    matchingHpaProteinExpression <- shiny.huge.hpaProteinExpession[symbol == matchingGene$symbol]

    scaledGtexValues <- matchingGtexExpression[,list(tissue = tissue, value = tpm_scaled)]

    scaledHpaRnaValues <- matchingHpaRnaExpression[,list(tissue = tissue, value = tpm_scaled)]

    hpaProteinValues <- matchingHpaProteinExpression[,list(value = max(level)),by = tissue]

    output$modalGTExScaledExpression <- modalExpressionPlot(scaledGtexValues, expressionFilter, paste(selectedSymbol, "GTEx data (scaled TPM)", sep = ": "))
    output$modalHpaRnaScaledExpression <- modalExpressionPlot(scaledHpaRnaValues, expressionFilter, paste(selectedSymbol, "HPA RNA data (scaled TPM)", sep = ":"))
    output$modalHpaProteinExpression <- modalExpressionPlot(hpaProteinValues, expressionFilter, paste(selectedSymbol, "HPA Protein data (levels)", sep = ":"))
    output$modalMutationTypes <- renderTable(uniqueMutations, spacing = "xs")

  })

}