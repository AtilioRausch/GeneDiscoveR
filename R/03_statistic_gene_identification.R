gene_identification_by_phenotype <- function(formula, PhenoRobject = NULL, statistic = "Fisher", name = "PerType", cores = 1) {
    if (!is.null(PhenoRobject$Identification)) {
        if (any(PhenoRobject$Identification$name == name & PhenoRobject$Identification$statistic == statistic)) {
            stop("Error: An execution with the same name and statistic already exists in PhenoRobject$Identification.")
        }
    }
    if (statistic == "Fisher") {
        # Comprueba si ambos elementos de la fórmula están en PhenoRobject$Phenotypes
        formula <- as.character(formula)
        formula_elements <- str_split(formula, pattern = "~")
        formula_elements <- c(formula[[2]], formula[[3]])
        # Check if any element of the formula contains "+"
        if (str_detect(formula_elements[1], pattern = "\\+")) {
            selected_elements <- str_split(formula_elements[1], "\\+", simplify = TRUE)
            selected_elements <- str_remove_all(selected_elements, pattern = " ")
            response <- c()
            for (i in 1:length(selected_elements)) {
                if (selected_elements[i] %in% names(PhenoRobject$Phenotypes)) {
                    response <- union(response, unlist(PhenoRobject$Phenotypes[selected_elements[i]]))
                } else {
                    stop("Error: '", selected_elements[i], "' is not a valid column name in PhenoRobject$Phenotypes.")
                }
            }
        } else {
            if (formula_elements[1] %in% names(PhenoRobject$Phenotypes)) {
                response <- as.vector(unlist(PhenoRobject$Phenotypes[formula_elements[1]]))
            } else {
                stop("Error: '", formula_elements[1], "' is not a valid column name in PhenoRobject$Phenotypes.")
            }
        }
        if (str_detect(formula_elements[2], pattern = "\\+")) {
            selected_elements <- str_split(formula_elements[2], "\\+", simplify = TRUE)
            selected_elements <- str_remove_all(selected_elements, pattern = " ")
            predictor <- c()
            for (i in 1:length(selected_elements)) {
                if (selected_elements[i] %in% names(PhenoRobject$Phenotypes)) {
                    predictor <- union(predictor, unlist(PhenoRobject$Phenotypes[selected_elements[i]]))
                } else {
                    stop("Error: '", selected_elements[i], "' is not a valid column name in PhenoRobject$Phenotypes.")
                }
            }
        } else {
            if (formula_elements[2] %in% names(PhenoRobject$Phenotypes)) {
                predictor <- as.vector(unlist(PhenoRobject$Phenotypes[formula_elements[2]]))
            } else {
                stop("Error: '", formula_elements[2], "' is not a valid column name in PhenoRobject$Phenotypes.")
            }
        }

        nameColumn1 <- paste0("N_", str_replace_all(pattern = " ", replacement = "_", formula[[3]]))
        nameColumn2 <- paste0("N_", str_replace_all(pattern = " ", replacement = "_", formula[[2]]))
        # Realiza la prueba exacta de Fisher
        # PhenoRobject$RunActive$N0Active <- PhenoRobject$RunActive$N0Active %>%
        #     rowwise() %>%
        #     mutate(
        #         nameColumn1 = sum(!is.na(c_across(all_of(predictor)))),
        #         nameColumn2 = sum(!is.na(c_across(all_of(response))))
        #     )
        nameColumn3 <- paste0("fisherResult", name)
        if (cores == 1) {
            PhenoRobject$RunActive$N0Active <- process_rows(PhenoRobject$RunActive$N0Active, predictor, response, nameColumn1, nameColumn2)
            PhenoRobject$RunActive$N0Active <- process_rows2(PhenoRobject$RunActive$N0Active, predictor, response, nameColumn1, nameColumn2, nameColumn3)
        } else {
            is_windows <- .Platform$OS.type == "windows"
            numCores <- min(cores, parallel::detectCores())

            split_data <- split(PhenoRobject$RunActive$N0Active, cut(seq(nrow(PhenoRobject$RunActive$N0Active)), numCores, labels = FALSE))
            # Procesar cada parte en paralelo
            tmpData <- mclapply(split_data, process_rows, predictor, response, nameColumn1, nameColumn2, mc.cores = numCores)
            PhenoRobject$RunActive$N0Active <- bind_rows(tmpData)
            rm(tmpData, split_data)

            split_data <- split(PhenoRobject$RunActive$N0Active, cut(seq(nrow(PhenoRobject$RunActive$N0Active)), numCores, labels = FALSE))
            tmpData <- mclapply(split_data, process_rows2, predictor, response, nameColumn1, nameColumn2, nameColumn3, mc.cores = numCores)
            PhenoRobject$RunActive$N0Active <- bind_rows(tmpData)
            rm(tmpData, split_data)
        }
        PhenoRobject$RunActive$N0Active <- PhenoRobject$RunActive$N0Active %>%
            rowwise() %>%
            mutate(
                !!paste0("pvalueFisher", name) := (!!sym(nameColumn3))[["p.value"]],
                !!paste0("oddsRatioFisher", name) := (!!sym(nameColumn3))[["estimate"]]
            ) %>%
            rename(!!nameColumn1 := nameColumn1, !!nameColumn2 := nameColumn2)
    }
    if (is.null(PhenoRobject$Identification)) {
        PhenoRobject$Identification[[1]] <- PhenoRIdentification(
            name, statistic, formula,
            Phenotypes = list("response" = response, "predictor" = predictor),
            columns = c(nameColumn1, nameColumn2, nameColumn3, paste0("pvalueFisher", name), paste0("oddsRatioFisher", name))
        )
    } else {
        PhenoRobject$Identification[[length(PhenoRobject$Identification) + 1]] <- PhenoRIdentification(
            name, statistic, formula,
            Phenotypes = list("response" = response, "predictor" = predictor),
            columns = c(nameColumn1, nameColumn2, nameColumn3, paste0("pvalueFisher", name), paste0("oddsRatioFisher", name))
        )
    }
    return(PhenoRobject)
}

select_genes_by_phenotype <- function(PhenoRobject, name = NULL, pvalue = 0.05, oddsRatio = 1, sign = ">=") {
    if (is.null(PhenoRobject)) {
        stop("Error: 'PhenoRobject' cannot be NULL.")
    }
    if (is.null(name)) {
        stop("Error: 'name' cannot be NULL.")
    }
    if (!any(endsWith(names(PhenoRobject$RunActive$N0Active), name))) {
        stop("Error: None of the columns in N0Active ends with '", name, "'.")
    }
    if (sign == ">=") {
        filterTable <- PhenoRobject$RunActive$N0Active %>%
            filter(!!sym(paste0("pvalueFisher", name)) <= pvalue && !!sym(paste0("oddsRatioFisher", name)) >= oddsRatio)
    } else if (sign == "<=") {
        filterTable <- PhenoRobject$RunActive$N0Active %>%
            filter(!!sym(paste0("pvalueFisher", name)) <= pvalue && !!sym(paste0("oddsRatioFisher", name)) <= oddsRatio)
    } else if (sign == "<") {
        filterTable <- PhenoRobject$RunActive$N0Active %>%
            filter(!!sym(paste0("pvalueFisher", name)) <= pvalue && !!sym(paste0("oddsRatioFisher", name)) < oddsRatio)
    } else if (sign == ">") {
        filterTable <- PhenoRobject$RunActive$N0Active %>%
            filter(!!sym(paste0("pvalueFisher", name)) <= pvalue && !!sym(paste0("oddsRatioFisher", name)) > oddsRatio)
    } else {
        stop("Error: 'sign' must be '>=', '<=', '<' or '>'.")
    }
    if (is.null(PhenoRobject$FilteredGenes)) {
        PhenoRobject$FilteredGenes[[1]] <- PhenoRFilteredGenes(table = filterTable, name = name, pvalue = pvalue, oddsRatio = oddsRatio, sign = sign)
    } else {
        PhenoRobject$FilteredGenes[[length(PhenoRobject$FilteredGenes) + 1]] <- PhenoRFilteredGenes(table = filterTable, name = name, pvalue = pvalue, oddsRatio = oddsRatio, sign = sign)
    }
    return(PhenoRobject)
}
