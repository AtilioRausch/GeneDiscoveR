#' Server function for GeneDiscoveR
#'
#' This function defines the server logic for the GeneDiscoveR application.
#' It handles the rendering of the menu, the selection of HOGs in the volcano plot,
#' and the filtering of data based on selected genes.
#'
#' @param input The input values from the UI.
#' @param output The output objects to be rendered in the UI.
#' @param session The session object for the Shiny application.
#'
#' @return None
#' @keywords internal
#' @examples
#' \dontrun{
#' .server_genediscover(input, output, session)
#' }
#' @import shinydashboard
#' @importFrom DT renderDT datatable
#' @importFrom plotly renderPlotly ggplotly event_data
#' @export
.server_genediscover <- function(input, output, session) {
    options(warn = -1)
    output$menu <- renderMenu({
        sidebarMenu(
            menuItem("Plot",
                tabName = "plot", icon = icon("bar-chart"), selected = TRUE,
                startExpanded = TRUE
            ),
            menuItem("Table of OGs", tabName = "tableOG", icon = icon("th"))
        )
    })
    isolate({
        updateTabItems(session, "tabs", "m2")
    })

    .set_ggplot2_theme()
    observeEvent(event_data("plotly_selected", source = "volcano"), {
        event <- event_data("plotly_selected", source = "volcano")
        if (!is.null(event) && length(event$pointNumber) > 0) {
            showNotification("Indication: The HOGs have been selected, you can view them in the ‘Table of orthogroups’ tab.", type = "default", duration = 10)
        }
    })

    selectedData <- reactive({
        event <- event_data("plotly_selected", source = "volcano")
        if (!is.null(event) && length(event$pointNumber) > 0) {
            dataselect <- suppressMessages(read_tsv(file = paste0(tempdir(), "/data-select.tmp")))
            dataunselect <- suppressMessages(read_tsv(file = paste0(tempdir(), "/data-unselect.tmp")))
            excludecols <- c("original_index", "contains-gene", "log-odds-ratio")
            excludecols <- c(excludecols, grep("^fisherResult", names(dataselect), value = TRUE))
            selected_data <- lapply(1:nrow(event), function(i) {
                if (event$curveNumber[i] == 0) {
                    dataunselect[event$pointNumber[i] + 1, -which(names(dataunselect) %in% excludecols)]
                } else if (event$curveNumber[i] == 1) {
                    dataselect[event$pointNumber[i] + 1, -which(names(dataselect) %in% excludecols)]
                }
            })
            # Convertir la lista a un dataframe
            do.call(rbind, selected_data)
        } else {
            data.frame()
        }
    })

    # Renderiza la tabla con los datos seleccionados
    output$table <- renderDT({
        datatable(
            selectedData(),
            options = list(
                scrollX = TRUE,
                pageLength = -1,
                scrollY = "calc(100vh - 200px)", # Ajusta según el espacio ocupado por otros elementos
                buttons = c("copy", "csv", "excel"),
                dom = "Bfrtip"
            ),
            filter = "top",
            extensions = "Buttons"
        )
    })
    # Volcano plot event
    observeEvent(input$submit, {
        # Verifica si el dataframe existe en el entorno global
        if (exists(input$dataframeName)) {
            data <- get(input$dataframeName)
            name <- input$name
            # Obtiene el dataframe por su nombre
            if (class(data) == "GeneDiscoveR" && is.character(name) && name %in% get_names_identification(GeneDiscoveRobject)) {
                showNotification("Success: The GeneDiscoveR object and the identification name were entered correctly. Please wait a moment!.", type = "message", duration = 10)

                GeneDiscoveRidentification <- .get_identification(GeneDiscoveRobject = GeneDiscoveRobject, name = name)
                data <- data$RunActive$N0Active

                saveRDS(GeneDiscoveRidentification, file = paste0(tempdir(), "/data.rds"))

                data <- data %>%
                    mutate(
                        "log-odds-ratio" = case_when(
                            !!sym(GeneDiscoveRidentification$columns[5]) == 0 ~ -5,
                            !!sym(GeneDiscoveRidentification$columns[5]) == Inf ~ 5,
                            is.finite(!!sym(GeneDiscoveRidentification$columns[5])) ~ log(!!sym(GeneDiscoveRidentification$columns[5]))
                        )
                    )
                data <- data %>%
                    rename(
                        "p-value" = GeneDiscoveRidentification$columns[4],
                        "odds-ratio" = GeneDiscoveRidentification$columns[5]
                    ) %>%
                    mutate("contains-gene" = FALSE)
                data$original_index <- seq_len(nrow(data))
                suppressMessages(write_tsv(data[data$`contains-gene` == TRUE, ], file = paste0(tempdir(), "/data-select.tmp")))
                suppressMessages(write_tsv(data[data$`contains-gene` == FALSE, ], file = paste0(tempdir(), "/data-unselect.tmp")))
                suppressMessages(write_tsv(data, file = paste0(tempdir(), "/data.tmp")))
                data <- data[order(data$original_index), ]
                output$plot <- renderPlotly({
                    g1 <- ggplot(data) +
                        geom_point(aes(x = -log(`p-value`), y = `log-odds-ratio`)) +
                        coord_flip() +
                        geom_hline(yintercept = 0, linetype = "dotted", col = "black") +
                        geom_vline(xintercept = 2.9957, linetype = "dotted", col = "black") +
                        annotate("text", x = 3.1, y = 1, label = "p-value <= 0.05", vjust = 1.5, color = "black") +
                        annotate("text", x = 8, y = 1, label = "OddRatio >= 1", color = "black") +
                        ylab("log(Odds Ratio)") +
                        xlab("-log(p-value)") +
                        scale_color_jama()
                    p <- ggplotly(g1, source = "volcano", dynamicTicks = TRUE) %>% toWebGL()
                })
                showNotification("Indication: you can select the HOGs in the Volcano plot!", type = "message", duration = 10)
            } else {
                if (class(data) == "GeneDiscoveR") {
                    showNotification("Error: The entered name is not the name of a GeneDiscoveR object!", type = "error")
                } else if (is.character(name) && name %in% get_names_identification(GeneDiscoveRobject)) {
                    showNotification("Error: The entered name does not correspond to an identification!", type = "error")
                }
            }
        }
    })

    observeEvent(input$submitgenes, {
        if (is.character(input$genes) && input$genes != "") {
            data <- suppressMessages(read_tsv(file = paste0(tempdir(), "/data.tmp")))
            GeneDiscoveRidentification <- readRDS(paste0(tempdir(), "/data.rds"))
            splitGenes <- str_split(input$genes, ",")[[1]]
            splitGenes <- trimws(splitGenes)

            showNotification("Success: The GeneDiscoveR object and the identification name were entered correctly. Please wait a moment!.", type = "message", duration = 10)

            data <- data %>%
            mutate(
                "contains-gene" = rowSums(sapply(splitGenes, function(gene) apply(data, 1, function(row) any(grepl(gene, row))))) > 0
            )
            
            data$original_index <- seq_len(nrow(data))
            suppressMessages(write_tsv(data[data$`contains-gene` == TRUE, ], file = paste0(tempdir(), "/data-select.tmp")))
            suppressMessages(write_tsv(data[data$`contains-gene` == FALSE, ], file = paste0(tempdir(), "/data-unselect.tmp")))

            if (any(data$`contains-gene` == TRUE)) {
            output$plot <- renderPlotly({
                g1 <- ggplot(data = data) +
                geom_point(aes(x = -log(`p-value`), y = `log-odds-ratio`, color = factor(`contains-gene`, labels = c("No", "Yes")))) +
                coord_flip() +
                geom_hline(yintercept = 0, linetype = "dotted", col = "black") +
                geom_vline(xintercept = 2.9957, linetype = "dotted", col = "black") +
                annotate("text", x = 3.1, y = 1, label = "p-value <= 0.05", vjust = 1.5, color = "black") +
                annotate("text", x = 8, y = 1, label = "Odds Ratio >= 1", color = "black") +
                ylab("log(Odds Ratio)") +
                xlab("-log(p-value)") +
                scale_color_jama() +
                guides(color = guide_legend(title = "Selected"))
                p <- ggplotly(g1, source = "volcano", dynamicTicks = TRUE) %>% toWebGL()
            })
            }
        }
})
}

