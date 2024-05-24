library(shiny)
library(plotly)
library(DT)
library(shinydashboard)

launch_web_app <- function() {
  ui <- dashboardPage(
    skin = "black",
    dashboardHeader(title = span(img(src = "/home/atilio/Escritorio/Git/GeneDiscoveR/R/logo.png", height = 35), "GeneDiscoveR")),
    dashboardSidebar(
      sidebarMenu(id = "tabs"),
      sidebarMenuOutput("menu")
    ),
    dashboardBody(
      tabItems(
        # First tab content
        tabItem(
          tabName = "plot",
          fluidRow(
            column(
              width = 6,
              box(
                title = "Instructions", status = "info", solidHeader = TRUE, width = NULL
              ),
              box(
                title = "Step 1: inputs", status = "warning", solidHeader = TRUE, width = NULL,
                textInput("dataframeName", "Enter the name of a GeneDiscoveR object:", value = "GeneDiscoveRobject"),
                textInput("name", "Enter the name of an identification:", value = "Self-compatible"),
                actionButton("submit", "Load", class = "btn-primary", icon = icon("play")),
              ),
              box(
                title = "Step 3: gene selector", status = "primary", solidHeader = TRUE, width = NULL,
                textInput("genes", 'Enter gene IDs separated by ",":', value = "Mp1g25500.1"),
                actionButton("submitgenes", "Load genes", class = "btn-primary", icon = icon("fa-dna")),
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
  shinyApp(ui, server)
}
