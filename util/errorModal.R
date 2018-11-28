util.showErrorModal <- function (errorMessage, session) {

  showModal(modalDialog(
    tags$div(errorMessage),
    title = "ERROR",
    easyClose = TRUE
  ))

}
