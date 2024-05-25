#' Launches the GeneDiscoveR Web App
#'
#' This function creates and launches the GeneDiscoveR web application. The web app provides a user interface for analyzing gene data using the GeneDiscoveR package.
#'
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' .launch_genediscover_web_app()
#' }
#' @import shiny shinydashboard plotly DT
.launch_genediscover_web_app <- function() {
  ui <- dashboardPage(
    skin = "black",
    dashboardHeader(
      # Título con enlace
      title = "GeneDiscoveR Web App"
    ),
    dashboardSidebar(
      sidebarMenu(id = "tabs"),
      sidebarMenuOutput("menu")
    ),
    dashboardBody(
      tags$style(HTML("
        body, input, select, textarea, .skin-black, .content-wrapper, .right-side {
          font-size: 0.75vw  /* Tamaño de la fuente relativo al ancho de la pantalla */
        }
        .box-title, h1, h2, h3, h4, h5, h6, .sidebar-menu, .main-sidebar{
          font-size: 1vw; /* Tamaño de los títulos relativo al ancho de la pantalla */
        }
      ")),
      tabItems(
        # First tab content
        tabItem(
          tabName = "plot",
          fluidRow(
            column(
              width = 6,
              box(
                title = "Instructions", status = "info", solidHeader = TRUE, width = NULL,
                # Agrega el logo a la derecha
                tags$div(
                  tags$img(src = "www/nombre_de_tu_imagen.png", height = "auto", width = "30%", style = "float: right; padding-left: 20px;"),
                  # Agrega la lista de instrucciones a la izquierda
                  tags$p("Welcome to the GeneDiscoveR work dashboard. To analyze your data, please follow these instructions:"),
                  tags$ul(
                    tags$li("Step 1: Enter the name of the GeneDiscoveR object and the name of the identification made with it in the “Step 1: Inputs” window. Then click Load."),
                    tags$li("Step 2: The associated Volcano plot will be loaded in “Step 2: Volcano plot”. In the plot, you can select HOGs using the Box select option and view them in the Table of HOGs section."),
                    tags$li("Step 3: Specify the names of the genes to be sectioned. For multiple selections, separate the IDs with a comma (,). When you click Load genes, the Volcano plot will be updated to allow their identification."),
                    tags$p("At any time, after selecting with Box select, you can go to the Table of HOGs section to view the data table, filter it, and/or export it.")
                  )
                ),
                style = "overflow:auto;"
              ),
              box(
                title = "Step 1: Inputs", status = "warning", solidHeader = TRUE, width = NULL,
                textInput("dataframeName", "Enter the name of a GeneDiscoveR object:", value = "GeneDiscoveRobject"),
                textInput("name", "Enter the name of an identification:", value = "Self-compatible"),
                actionButton("submit", "Load", class = "btn-primary", icon = icon("play")),
              ),
              box(
                title = "Step 3: Gene selector", status = "primary", solidHeader = TRUE, width = NULL,
                textInput("genes", 'Enter gene IDs separated by ",":', value = "Mp1g25500.1"),
                actionButton("submitgenes", "Load genes", class = "btn-primary", icon = icon("dna")),
              ),
            ),
            box(
              title = "Step 2: Volcano plot", status = "primary", solidHeader = TRUE,
              collapsible = TRUE,
              plotlyOutput("plot", height = "80vh") # Set height to "auto" for adjustable height
            ),
          ),
        ),

        # Second tab content
        tabItem(
          tabName = "tableHOG",
          dataTableOutput("table")
        )
      )
    )
  )
  shinyApp(ui, .server_genediscover)
}
