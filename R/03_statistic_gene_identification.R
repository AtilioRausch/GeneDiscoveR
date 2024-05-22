gene_identification_by_phenotype <- function(formula, GeneDiscoveRobject = NULL, statistic = "Fisher", name = "PerType", cores = 1) {
    if (!is.null(GeneDiscoveRobject$Identification)) {
        if (any(GeneDiscoveRobject$Identification$name == name & GeneDiscoveRobject$Identification$statistic == statistic)) {
            stop("Error: An execution with the same name and statistic already exists in GeneDiscoveRobject$Identification.")
        }
    }
    if (statistic == "Fisher") {
        # Comprueba si ambos elementos de la fórmula están en GeneDiscoveRobject$Phenotypes
        formula <- as.character(formula)
        formula_elements <- str_split(formula, pattern = "~")
        formula_elements <- c(formula[[2]], formula[[3]])
        # Check if any element of the formula contains "+"
        if (str_detect(formula_elements[1], pattern = "\\+")) {
            selected_elements <- str_split(formula_elements[1], "\\+", simplify = TRUE)
            selected_elements <- str_remove_all(selected_elements, pattern = " ")
            response <- c()
            for (i in 1:length(selected_elements)) {
                if (selected_elements[i] %in% names(GeneDiscoveRobject$Phenotypes)) {
                    response <- union(response, unlist(GeneDiscoveRobject$Phenotypes[selected_elements[i]]))
                } else {
                    stop("Error: '", selected_elements[i], "' is not a valid column name in GeneDiscoveRobject$Phenotypes.")
                }
            }
        } else {
            if (formula_elements[1] %in% names(GeneDiscoveRobject$Phenotypes)) {
                response <- as.vector(unlist(GeneDiscoveRobject$Phenotypes[formula_elements[1]]))
            } else {
                stop("Error: '", formula_elements[1], "' is not a valid column name in GeneDiscoveRobject$Phenotypes.")
            }
        }
        if (str_detect(formula_elements[2], pattern = "\\+")) {
            selected_elements <- str_split(formula_elements[2], "\\+", simplify = TRUE)
            selected_elements <- str_remove_all(selected_elements, pattern = " ")
            predictor <- c()
            for (i in 1:length(selected_elements)) {
                if (selected_elements[i] %in% names(GeneDiscoveRobject$Phenotypes)) {
                    predictor <- union(predictor, unlist(GeneDiscoveRobject$Phenotypes[selected_elements[i]]))
                } else {
                    stop("Error: '", selected_elements[i], "' is not a valid column name in GeneDiscoveRobject$Phenotypes.")
                }
            }
        } else {
            if (formula_elements[2] %in% names(GeneDiscoveRobject$Phenotypes)) {
                predictor <- as.vector(unlist(GeneDiscoveRobject$Phenotypes[formula_elements[2]]))
            } else {
                stop("Error: '", formula_elements[2], "' is not a valid column name in GeneDiscoveRobject$Phenotypes.")
            }
        }

        nameColumn1 <- paste0("N_", str_replace_all(pattern = " ", replacement = "_", formula[[3]]))
        nameColumn2 <- paste0("N_", str_replace_all(pattern = " ", replacement = "_", formula[[2]]))
        # Realiza la prueba exacta de Fisher
        # GeneDiscoveRobject$RunActive$N0Active <- GeneDiscoveRobject$RunActive$N0Active %>%
        #     rowwise() %>%
        #     mutate(
        #         nameColumn1 = sum(!is.na(c_across(all_of(predictor)))),
        #         nameColumn2 = sum(!is.na(c_across(all_of(response))))
        #     )
        nameColumn3 <- paste0("fisherResult", name)
        if (cores == 1) {
            GeneDiscoveRobject$RunActive$N0Active <- process_rows(GeneDiscoveRobject$RunActive$N0Active, predictor, response, nameColumn1, nameColumn2)
            GeneDiscoveRobject$RunActive$N0Active <- process_rows2(GeneDiscoveRobject$RunActive$N0Active, predictor, response, nameColumn1, nameColumn2, nameColumn3)
        } else {
            is_windows <- .Platform$OS.type == "windows"
            numCores <- min(cores, parallel::detectCores())

            split_data <- split(GeneDiscoveRobject$RunActive$N0Active, cut(seq(nrow(GeneDiscoveRobject$RunActive$N0Active)), numCores, labels = FALSE))
            # Procesar cada parte en paralelo
            tmpData <- mclapply(split_data, process_rows, predictor, response, nameColumn1, nameColumn2, mc.cores = numCores)
            GeneDiscoveRobject$RunActive$N0Active <- bind_rows(tmpData)
            rm(tmpData, split_data)

            split_data <- split(GeneDiscoveRobject$RunActive$N0Active, cut(seq(nrow(GeneDiscoveRobject$RunActive$N0Active)), numCores, labels = FALSE))
            tmpData <- mclapply(split_data, process_rows2, predictor, response, nameColumn1, nameColumn2, nameColumn3, mc.cores = numCores)
            GeneDiscoveRobject$RunActive$N0Active <- bind_rows(tmpData)
            rm(tmpData, split_data)
        }
        GeneDiscoveRobject$RunActive$N0Active <- GeneDiscoveRobject$RunActive$N0Active %>%
            rowwise() %>%
            mutate(
                !!paste0("pvalueFisher", name) := (!!sym(nameColumn3))[["p.value"]],
                !!paste0("oddsRatioFisher", name) := (!!sym(nameColumn3))[["estimate"]]
            ) %>%
            rename(!!nameColumn1 := nameColumn1, !!nameColumn2 := nameColumn2)
    }
    if (is.null(GeneDiscoveRobject$Identification)) {
        GeneDiscoveRobject$Identification[[1]] <- GeneDiscoveRIdentification(
            name, statistic, formula,
            Phenotypes = list("response" = response, "predictor" = predictor),
            columns = c(nameColumn1, nameColumn2, nameColumn3, paste0("pvalueFisher", name), paste0("oddsRatioFisher", name))
        )
    } else {
        GeneDiscoveRobject$Identification[[length(GeneDiscoveRobject$Identification) + 1]] <- GeneDiscoveRIdentification(
            name, statistic, formula,
            Phenotypes = list("response" = response, "predictor" = predictor),
            columns = c(nameColumn1, nameColumn2, nameColumn3, paste0("pvalueFisher", name), paste0("oddsRatioFisher", name))
        )
    }
    return(GeneDiscoveRobject)
}

select_genes_by_phenotype <- function(GeneDiscoveRobject, name = NULL, pvalue = 0.05, oddsRatio = 1, sign = ">=") {
    if (is.null(GeneDiscoveRobject)) {
        stop("Error: 'GeneDiscoveRobject' cannot be NULL.")
    }
    if (is.null(name)) {
        stop("Error: 'name' cannot be NULL.")
    }
    if (!any(endsWith(names(GeneDiscoveRobject$RunActive$N0Active), name))) {
        stop("Error: None of the columns in N0Active ends with '", name, "'.")
    }
    if (sign == ">=") {
        filterTable <- GeneDiscoveRobject$RunActive$N0Active %>%
            filter(!!sym(paste0("pvalueFisher", name)) <= pvalue && !!sym(paste0("oddsRatioFisher", name)) >= oddsRatio)
    } else if (sign == "<=") {
        filterTable <- GeneDiscoveRobject$RunActive$N0Active %>%
            filter(!!sym(paste0("pvalueFisher", name)) <= pvalue && !!sym(paste0("oddsRatioFisher", name)) <= oddsRatio)
    } else if (sign == "<") {
        filterTable <- GeneDiscoveRobject$RunActive$N0Active %>%
            filter(!!sym(paste0("pvalueFisher", name)) <= pvalue && !!sym(paste0("oddsRatioFisher", name)) < oddsRatio)
    } else if (sign == ">") {
        filterTable <- GeneDiscoveRobject$RunActive$N0Active %>%
            filter(!!sym(paste0("pvalueFisher", name)) <= pvalue && !!sym(paste0("oddsRatioFisher", name)) > oddsRatio)
    } else {
        stop("Error: 'sign' must be '>=', '<=', '<' or '>'.")
    }
    if (is.null(GeneDiscoveRobject$FilteredGenes)) {
        GeneDiscoveRobject$FilteredGenes[[1]] <- GeneDiscoveRFilteredGenes(table = filterTable, name = name, pvalue = pvalue, oddsRatio = oddsRatio, sign = sign)
    } else {
        GeneDiscoveRobject$FilteredGenes[[length(GeneDiscoveRobject$FilteredGenes) + 1]] <- GeneDiscoveRFilteredGenes(table = filterTable, name = name, pvalue = pvalue, oddsRatio = oddsRatio, sign = sign)
    }
    return(GeneDiscoveRobject)
}
