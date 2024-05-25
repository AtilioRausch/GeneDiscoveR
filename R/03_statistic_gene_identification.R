#' Gene Identification by Phenotype
#'
#' This function performs gene identification based on a given formula, using the Fisher statistic.
#' The identification permits the comparison of the relationship between the response and predictor variables.
#' This function is used to separate the orthogroups by their enrichment in species with analized phenotypes.
#' In the case of multiple response or predictor variables, the formula must be separated by a plus sign (+).
#' For example, if the formula is "one_in_specialized_cell ~ many_in_all_cells", the function will apply the Fisher test to distinguish the orthogroups that are present in the species with the phenotype "one_in_specialized_cell" and absent in the species with the phenotype "many_in_all_cells".
#' and the orthogroups that are present in the species with the phenotype "many_in_all_cells" and absent in the species with the phenotype "one_in_specialized_cell".
#' The function needs the previous selection of species by phenotype with \code{\link{select_species_by_phenotype}}.
#' @param formula The formula specifying the relationship between the response and predictor variables.
#' @param GeneDiscoveRobject An optional GeneDiscoveR object to store the identification results.
#' @param statistic The statistical test to be used for gene identification. Default is "Fisher".
#' @param name The name of the identification execution.
#' @param cores The number of cores to be used for parallel processing. Default is 1.
#'
#' @return The updated GeneDiscoveR object with the identification results.
#'
#' @examples
#'
#' # Create a GeneDiscoveR object
#' N0sDir <- system.file("extdata", "N0-1dot3-6", package = "GeneDiscoveR")
#' overallsDir <- system.file("extdata", "Comparatives-1dot3-6", package = "GeneDiscoveR")
#' dataFile <- system.file("extdata", "annotatedCDSs.tsv", package = "GeneDiscoveR")
#' minInflation <- 1.3
#' maxInflation <- 6
#' stepInflation <- 0.1
#'
#' GeneDiscoveRobject <- GeneDiscoveR(overallsDir = overallsDir, N0sDir = N0sDir, dataFile = dataFile, minInflation = minInflation, maxInflation = maxInflation, stepInflation = stepInflation)
#'
#' # Set active run
#' GeneDiscoveRobject <- set_run_active(GeneDiscoveRobject, InflationValue = 1.8)
#'
#' # Select species by phenotype
#' GeneDiscoveRobject <- select_species_by_phenotype(GeneDiscoveRobject = GeneDiscoveRobject, columnPhenotype = "Oil-body-type", columnID = "OrthofinderID", type = "one_in_specialized_cell")
#' GeneDiscoveRobject <- select_species_by_phenotype(GeneDiscoveRobject = GeneDiscoveRobject, columnPhenotype = "Oil-body-type", columnID = "OrthofinderID", type = "many_in_all_cells")
#' GeneDiscoveRobject <- select_species_by_phenotype(GeneDiscoveRobject = GeneDiscoveRobject, columnPhenotype = "Oil-body-type", columnID = "OrthofinderID", type = "noneOB")
#'
#' # Gene identification by phenotype
#' # Identifies genes that are present in the species with the phenotype "one_in_specialized_cell" and absent in the species with the phenotype "many_in_all_cells".
#' GeneDiscoveRobject <- gene_identification_by_phenotype(GeneDiscoveRobject = GeneDiscoveRobject, formula = as.formula("one_in_specialized_cell ~ many_in_all_cells"), statistic = "Fisher", name = "PerOBtype", cores = 8)
#' # Identifies genes that are present in the species with the phenotype "one_in_specialized_cell + many_in_all_cells" and absent in the species with the phenotype "noneOB".
#' GeneDiscoveRobject <- gene_identification_by_phenotype(GeneDiscoveRobject = GeneDiscoveRobject, formula = as.formula("noneOB ~ one_in_specialized_cell + many_in_all_cells"), statistic = "Fisher", name = "OBpresence", cores = 8)
#'
#' @export
#' @import dplyr stringr
#' @import doParallel
#' @import parallel
gene_identification_by_phenotype <- function(formula, GeneDiscoveRobject = NULL, statistic = "Fisher", name = "PerType", cores = 1) {
    if (!is.null(GeneDiscoveRobject$Identification)) {
        if (any(GeneDiscoveRobject$Identification$name == name & GeneDiscoveRobject$Identification$statistic == statistic)) {
            stop("Error: An execution with the same name and statistic already exists in GeneDiscoveRobject$Identification.")
        }
    }
    if (statistic == "Fisher") {
        formula <- as.character(formula)
        formula_elements <- str_split(formula, pattern = "~")
        formula_elements <- c(formula[[2]], formula[[3]])

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

        nameColumn3 <- paste0("fisherResult", name)
        if (cores == 1) {
            GeneDiscoveRobject$RunActive$N0Active <- .sum_per_phenotype(GeneDiscoveRobject$RunActive$N0Active, predictor, response, nameColumn1, nameColumn2)
            GeneDiscoveRobject$RunActive$N0Active <- .fisher_per_row(GeneDiscoveRobject$RunActive$N0Active, predictor, response, nameColumn1, nameColumn2, nameColumn3)
        } else {
            is_windows <- .Platform$OS.type == "windows"
            numCores <- min(cores, detectCores())

            split_data <- split(GeneDiscoveRobject$RunActive$N0Active, cut(seq(nrow(GeneDiscoveRobject$RunActive$N0Active)), numCores, labels = FALSE))

            tmpData <- mclapply(split_data, .sum_per_phenotype, predictor, response, nameColumn1, nameColumn2, mc.cores = numCores)
            GeneDiscoveRobject$RunActive$N0Active <- bind_rows(tmpData)
            rm(tmpData, split_data)

            split_data <- split(GeneDiscoveRobject$RunActive$N0Active, cut(seq(nrow(GeneDiscoveRobject$RunActive$N0Active)), numCores, labels = FALSE))
            tmpData <- mclapply(split_data, .fisher_per_row, predictor, response, nameColumn1, nameColumn2, nameColumn3, mc.cores = numCores)
            GeneDiscoveRobject$RunActive$N0Active <- bind_rows(tmpData)
            rm(tmpData, split_data)
        }
        GeneDiscoveRobject$RunActive$N0Active <- GeneDiscoveRobject$RunActive$N0Active %>%
            rowwise() %>%
            mutate(
                !!paste0("pvalueFisher", name) := (!!sym(nameColumn3))[["p.value"]],
                !!paste0("oddsRatioFisher", name) := (!!sym(nameColumn3))[["estimate"]]
            )
        if (!any("nameColumn1" %in% colnames(GeneDiscoveRobject$RunActive$N0Active))) {
            GeneDiscoveRobject$RunActive$N0Active <- GeneDiscoveRobject$RunActive$N0Active %>%
                rename(!!nameColumn1 := nameColumn1)
        }
        if (!any("nameColumn2" %in% colnames(GeneDiscoveRobject$RunActive$N0Active))) {
            GeneDiscoveRobject$RunActive$N0Active <- GeneDiscoveRobject$RunActive$N0Active %>%
                rename(!!nameColumn2 := nameColumn2)
        }
    }
    if (is.null(GeneDiscoveRobject$Identification)) {
        GeneDiscoveRobject$Identification[[1]] <- .GeneDiscoveRIdentification(
            name, statistic, formula,
            Phenotypes = list("response" = response, "predictor" = predictor),
            columns = c(nameColumn1, nameColumn2, nameColumn3, paste0("pvalueFisher", name), paste0("oddsRatioFisher", name))
        )
    } else {
        GeneDiscoveRobject$Identification[[length(GeneDiscoveRobject$Identification) + 1]] <- .GeneDiscoveRIdentification(
            name, statistic, formula,
            Phenotypes = list("response" = response, "predictor" = predictor),
            columns = c(nameColumn1, nameColumn2, nameColumn3, paste0("pvalueFisher", name), paste0("oddsRatioFisher", name))
        )
    }
    return(GeneDiscoveRobject)
}

#' Select genes by phenotype
#'
#' This function filter the previuos identification with \code{\link{gene_identification_by_phenotype}}.
#' Permit to filter genes by p-value and odds ratio for found orthogroups enriched in genes from species with one phenotype and not in other.
#' With the \code{formula("firstPhenotype ~ secondPhenotype")}.
#' If odds ratio is greater than 1, the orthogroups are enriched in genes from species with the first phenotype.
#' If odds ratio is less than 1, the orthogroups are enriched in genes from species with the second phenotype.
#'
#' @param GeneDiscoveRobject The GeneDiscoveR object containing the necessary data.
#' @param name The name of the identification execution.
#' @param pvalue The p-value threshold for filtering genes.
#' @param oddsRatio The odds ratio threshold for filtering genes.
#' @param sign The sign of the comparison operator for odds ratio filtering genes. Must be one of '>=', '<=', '<', or '>'.
#'
#' @return The updated GeneDiscoveRobject with the filtered genes added to the FilteredGenes list.
#'
#' @examples
#' # Create a GeneDiscoveR object
#' N0sDir <- system.file("extdata", "N0-1dot3-6", package = "GeneDiscoveR")
#' overallsDir <- system.file("extdata", "Comparatives-1dot3-6", package = "GeneDiscoveR")
#' dataFile <- system.file("extdata", "annotatedCDSs.tsv", package = "GeneDiscoveR")
#' minInflation <- 1.3
#' maxInflation <- 6
#' stepInflation <- 0.1
#'
#' GeneDiscoveRobject <- GeneDiscoveR(overallsDir = overallsDir, N0sDir = N0sDir, dataFile = dataFile, minInflation = minInflation, maxInflation = maxInflation, stepInflation = stepInflation)
#'
#' # Set active run
#' GeneDiscoveRobject <- set_run_active(GeneDiscoveRobject, InflationValue = 1.8)
#'
#' # Select species by phenotype
#' GeneDiscoveRobject <- select_species_by_phenotype(GeneDiscoveRobject = GeneDiscoveRobject, columnPhenotype = "Oil-body-type", columnID = "OrthofinderID", type = "one_in_specialized_cell")
#' GeneDiscoveRobject <- select_species_by_phenotype(GeneDiscoveRobject = GeneDiscoveRobject, columnPhenotype = "Oil-body-type", columnID = "OrthofinderID", type = "many_in_all_cells")
#'
#' # Gene identification by phenotype
#' GeneDiscoveRobject <- gene_identification_by_phenotype(GeneDiscoveRobject = GeneDiscoveRobject, formula = as.formula("one_in_specialized_cell ~ many_in_all_cells"), statistic = "Fisher", name = "PerOBtype", cores = 8)
#'
#' # Select genes by phenotype
#' GeneDiscoveRobject <- select_genes_by_phenotype(GeneDiscoveRobject, pvalue = 0.05, oddsRatio = 1, sign = ">=", name = "PerOBtype")
#' GeneDiscoveRobject <- select_genes_by_phenotype(GeneDiscoveRobject, pvalue = 0.05, oddsRatio = 1, sign = "<=", name = "PerOBtype")
#'
#' # View filtered genes table
#'
#' filteredGenes <- get_filtered_genes_table(GeneDiscoveRobject, name = "PerOBtype", pvalue = 0.05, oddsRatio = 1, sign = ">=")
#' # Output: A table with the filtered orthogroups enrichment in genes from species with phenotype "one_in_specialized_cell".
#' filteredGenes <- get_filtered_genes_table(GeneDiscoveRobject, name = "PerOBtype", pvalue = 0.05, oddsRatio = 1, sign = "<=")
#' # Output: A table with the filtered orthogroups enrichment in genes from species with phenotype "many_in_all_cells".
#'
#' @export
#' @import dplyr
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
        GeneDiscoveRobject$FilteredGenes[[1]] <- .GeneDiscoveRFilteredGenes(table = filterTable, name = name, pvalue = pvalue, oddsRatio = oddsRatio, sign = sign)
    } else {
        GeneDiscoveRobject$FilteredGenes[[length(GeneDiscoveRobject$FilteredGenes) + 1]] <- .GeneDiscoveRFilteredGenes(table = filterTable, name = name, pvalue = pvalue, oddsRatio = oddsRatio, sign = sign)
    }
    return(GeneDiscoveRobject)
}
