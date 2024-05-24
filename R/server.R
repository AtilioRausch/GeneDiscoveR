server <- function(input, output, session) {
    # Side bar
    output$menu <- renderMenu({
        sidebarMenu(
            menuItem("Plot",
                tabName = "plot", icon = icon("bar-chart"), selected = TRUE,
                startExpanded = TRUE
            ),
            menuItem("Table of HOGs selected", tabName = "tableHOG", icon = icon("th"))
        )
    })
    isolate({
        updateTabItems(session, "tabs", "m2")
    })

    observeEvent(event_data("plotly_selected", source = "volcano"), {
        event <- event_data("plotly_selected", source = "volcano")
        if (!is.null(event)) {
            showNotification("Indication: The HOGs have been selected, you can view them in the ‘Table of selected HOGs’ tab.", type = "default", duration = 5)
        }
    })

    selectedData <- reactive({
        event <- event_data("plotly_selected", source = "volcano")
        if (!is.null(event) && length(event$pointNumber) > 0) {
            selected_indices <- event$pointNumber + 1
            excluded_columns <- grepl("^fisherResult", names(data)) | grepl("^logodd", names(data))
            data[selected_indices, !excluded_columns]
        } else {
            data.frame()
        }
    })

    # Renderiza la tabla con los datos seleccionados
    output$table <- renderDataTable({
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
                showNotification("Success: The GeneDiscoveR object and the identification name were entered correctly. Please wait a moment!.", type = "message", duration = NULL)

                GeneDiscoveRidentification <- get_identification(GeneDiscoveRobject = GeneDiscoveRobject, name = name)
                data <- data$RunActive$N0Active

                write_tsv(data, file = paste0(tempdir(), "/data.txt"))
                saveRDS(GeneDiscoveRidentification, file = paste0(tempdir(), "/data.rds"))

                data <- data %>%
                    mutate(
                        logoddRatioFisher = case_when(
                            !!sym(GeneDiscoveRidentification$columns[5]) == 0 ~ -5,
                            !!sym(GeneDiscoveRidentification$columns[5]) == Inf ~ 5,
                            is.finite(!!sym(GeneDiscoveRidentification$columns[5])) ~ log(!!sym(GeneDiscoveRidentification$columns[5]))
                        )
                    )
                output$plot <- renderPlotly({
                    g1 <- ggplot() +
                        geom_point(aes(x = -log(data[[GeneDiscoveRidentification$columns[4]]]), y = data$logoddRatioFisher), color = "#696D7D") +
                        coord_flip() +
                        geom_hline(yintercept = 0, linetype = "dotted", col = "red") +
                        geom_vline(xintercept = 2.9957, linetype = "dotted", col = "red") +
                        annotate("text", x = 2.9957, y = 0, label = "p-value <= 0.05", vjust = 1.5, color = "black") +
                        annotate("text", x = 8, y = 1, label = "OddRatio >= 1", color = "black") +
                        ylab("") +
                        xlab("-log(p-value)")
                    ggplotly(g1, source = "volcano", dynamicTicks = TRUE)
                })
                showNotification("Indication: you can select the HOGs in the Volcano plot!", type = "message", duration = NULL)
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
            data <- read_tsv(file = paste0(tempdir(), "/data.txt"))
            GeneDiscoveRidentification <- readRDS(paste0(tempdir(), "/data.rds"))
            splitGenes <- str_split(input$genes, ",")[[1]]
            splitGenes <- trimws(splitGenes)

            showNotification("Success: The GeneDiscoveR object and the identification name were entered correctly. Please wait a moment!.", type = "message", duration = NULL)

            data <- data %>%
                mutate(
                    "contains-gene" = rowSums(sapply(splitGenes, function(gene) apply(data, 1, function(row) any(grepl(gene, row))))) > 0
                )
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
                )

            output$plot <- renderPlotly({
                g1 <- ggplot(data = data) +
                    geom_point(aes(x = -log(`p-value`), y = `log-odds-ratio`, color = `contains-gene`)) +
                    coord_flip() +
                    geom_hline(yintercept = 0, linetype = "dotted", col = "black") +
                    geom_vline(xintercept = 2.9957, linetype = "dotted", col = "black") +
                    annotate("text", x = 2.9957, y = 0, label = "p-value <= 0.05", vjust = 1.5, color = "black") +
                    annotate("text", x = 8, y = 1, label = "Odds Ratio >= 1", color = "black") +
                    ylab("log(Odds Ratio)") +
                    xlab("-log(p-value)")

                ggplotly(g1, source = "volcano", dynamicTicks = TRUE)
            })
        }
    })
}
