library(shiny)
library(plotly)
library(DT)
launch_web_app <- function() {
  # Definir la interfaz de usuario (UI)
  ui <- fluidPage(
    tags$head(
      tags$style(HTML("
       #myImage {
      display: block;
      margin-top: auto;
      margin-left: auto;
      margin-right: auto;
      max-width: 75%; /* Asegura que la imagen no sea más ancha que el contenedor */
      height: auto; /* Mantiene la relación de aspecto de la imagen */
    }
    .shiny-image-output {
        display: flex;
        justify-content: center;
        align-items: center;
        height: 100%;
      }
    ")),
      tags$script(HTML('
      $(document).on("shiny:connected", function() {
        var plotContainer = document.querySelector(".plot-container");
        var aspectRatio = 4 / 3;
        var height = plotContainer.clientHeight;
        var width = height * aspectRatio;
        plotContainer.style.width = width + "px";
      });
    ')),
    ),
    div(
      style = "display: flex; justify-content: space-around; align-items: center;",
      titlePanel("GeneDiscoveR: Interactive Volcano Plot and Table", windowTitle = "GeneDiscoveR: Interactive Volcano Plot and Table"),
      div(
        class = "container", # Estilo para centrar y ajustar la imagen
        uiOutput("myImage")
      ),
    ),
    div(
      style = "display: flex; justify-content: space-between;",
      div(
        style = "display: flex; justify-content: center; align-items: center; flex-wrap:wrap", # Estilo para maximizar el tamaño del gráfico
        textInput("dataframeName", "Enter the name of a GeneDiscoveR object:"),
        textInput("name", "Enter the name of an identification:"),
        actionButton("submit", "Load", class = "btn-primary"),
      )
    ),
    tabsetPanel(
      tabPanel(
        "Volcano Plot",
        div(
          class = "plot-container",
          style = "display: flex; justify-content: center; align-items: center; height: 100vh;", # Estilo para maximizar el tamaño del gráfico
          plotlyOutput("plot", height = "100%")
        )
      ),
      tabPanel(
        "Table of selected HOGs",
        dataTableOutput("table")
      )
    )
  )
  # fluidPage(
  #   #titlePanel("Gráfico Interactivo y Tabla"),
  #   fluidRow(
  #     column(
  #       width = 12, # Para pantallas grandes
  #       class = "col-md-6 col-sm-12", # 6 columnas en pantallas medianas y más grandes, 12 columnas en pantallas pequeñas
  #       plotlyOutput("plot", height = "100%")
  #     ),
  #     column(
  #       width = 12, # Para pantallas grandes
  #       class = "col-md-6 col-sm-12", # 6 columnas en pantallas medianas y más grandes, 12 columnas en pantallas pequeñas
  #       tableOutput("table")
  #     )
  #   )
  # )
  # Lanzar la aplicación Shiny
  shinyApp(ui = ui, server = server)
}
