util.copyToClipboardButton <- function (text, id) {
  return(renderUI({
    rclipButton(id, "Copy to clipboard", text, icon("clipboard"))
  }))
}
