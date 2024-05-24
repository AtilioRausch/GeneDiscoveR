library(shiny)
library(plotly)
library(DT)
library(shinydashboard)

launch_genediscover_web_app <- function() {
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
                  tags$img(src = "www/logo.png", height = "35", style = "float: right;"),
                  # Agrega la lista de instrucciones a la izquierda
                  tags$ul(
                    tags$li("Step 1: Enter the name of a GeneDiscoveR object and an identification."),
                    tags$li("Step 2: Click 'Load' to load the data."),
                    tags$li("Step 3: Enter gene IDs separated by ','."),
                    tags$li("Step 4: Click 'Load genes' to load the genes."),
                    tags$li("Step 5: View the Volcano plot on the right.")
                  )
                )
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
  shinyApp(ui, server_genediscover)
}
