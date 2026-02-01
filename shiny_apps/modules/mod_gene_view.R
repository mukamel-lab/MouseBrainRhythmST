# mod_gene_view.R

mod_gene_view_ui <- function(id) {
  ns <- NS(id)
  tagList(
    textInput(ns("gene"), "Gene"),
    verbatimTextOutput(ns("out"))
  )
}

mod_gene_view_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$out <- renderPrint(input$gene)
  })
}
