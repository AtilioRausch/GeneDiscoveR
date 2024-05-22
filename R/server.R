server <- function(input, output) {
    # Define the desired width
    # Crear el gráfico interactivo

    output$myImage <- renderUI({
        image <- reactive({
            outfile <- "/home/atilio/Escritorio/LiverwortGitHub/GeneDiscoveR/R/logo.png"
            list(
                src = outfile,
                contentType = "image/png",
                alt = "This is alternate text",
                id = "myImage"
            )
        })
        renderImage(
            {
                image()
            },
            deleteFile = FALSE
        )
    })


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

                # Crear la tabla que refleje la selección en el gráfico
                # output$table <- renderTable({
                #     event <- event_data("plotly_selected", source = "volcano") # Usar el mismo source que el gráfico interactivo anterior
                #     if (!is.null(event)) {
                #         selected_indices <- event$pointNumber + 1
                #         if (length(selected_indices) > 0) {
                #             selected_data <- data[selected_indices, c("OG", "Gene Tree Parent Clade", "MpTAKv6-Marchantia_polymorpha_rudelaris", "oddsRatioFisherOBpresence", "pvalueFisherOBpresence")]
                #             selected_data
                #         } else {
                #             NULL
                #         }
                #     } else {
                #         NULL
                #     }
                # })
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
            } else {
                if (class(data) == "GeneDiscoveR") {
                    showNotification("Error: The entered name is not the name of a GeneDiscoveR object!", type = "error")
                } else if (is.character(name) && name %in% get_names_identification(GeneDiscoveRobject)) {
                    showNotification("Error: The entered name does not correspond to an identification!", type = "error")
                }
                # Si el dataframe no existe, muestra una notificación de error
            }
        }
    })
}
