#' Convert a table to character vector
#'
#' This function takes a table as input and converts it to a character vector.
#'
#' @param table The table to be converted. It should be a data frame or a tibble.
#' @return A character vector representing the table.
#' @examples
#' table <- data.frame(x = c(1, 2, 3), y = c("a", "b", "c"))
#' .table_as_character(table)
#' # Output: [1] "1" "2" "3" "a" "b" "c"
#' @keywords internal
#' @import dplyr
#' @export
.table_as_character <- function(table = NULL) {
    if (is.null(table)) {
        stop("Error: 'table' cannot be NULL.")
    }

    result <- table %>%
        pull(.) %>%
        as.character()

    return(result)
}

#' Create a GeneDiscoveR object
#'
#' This function creates a GeneDiscoveR object with the specified parameters.
#' The data can be originated from OrthoFinder or InParanoiDB. Can be indicate the origin of the orthologs definition with the parameter \code{orthologsTool}.
#'
#' For OrthoFinder, the parameters \code{overallsDir}, \code{N0sDir}, \code{dataFile}, \code{uniqueInflation}, \code{minInflation}, \code{maxInflation}, and \code{stepInflation} are required.
#' The parameters of Inflation depend on the type of analysis to be performed with unique inflation or a range of inflations.
#'
#' For InParanoiDB, the parameters \code{pairSpeciesDir}, \code{principalSpecie}, and \code{principalSpeciePrefix} are required.
#' The comparison is made between the principal species and the other species. The principal species is the species that will be used as a reference.
#' The principal species name is \code{principalSpecie} and the principal species gene prefix is \code{principalSpeciePrefix}.
#'
#' @param overallsDir The directory path for the overalls (only for OrthoFinder run).
#' @param N0sDir The directory path for the N0s (only for OrthoFinder run).
#' @param dataFile The file path for the table with at least two columns (speciesID and n columns one per trait or phenotype) data.
#' @param annotationFile The file path for the annotation data.
#' @param uniqueInflation The unique inflation value to be used (only for OrthoFinder run).
#' @param minInflation The minimum inflation value for the range of inflations (only for OrthoFinder run).
#' @param maxInflation The maximum inflation value for the range of inflations (only for OrthoFinder run).
#' @param stepInflation The step size for the range of inflations (only for OrthoFinder run).
#' @param pairSpeciesDir The directory path for the pair species (only for InParanoiDB run).
#' @param principalSpecie The principal species (only for InParanoiDB run).
#' @param principalSpeciePrefix The principal species prefix (only for InParanoiDB run).
#' @param orthologsTool The tool used for orthologs identification, default is OrthoFinder (OrthoFinder or InParanoiDB).
#'
#' @return A GeneDiscoveR object with the specified parameters.
#' @export
#'
#' @examples
#' # Create a GeneDiscoveR for an execution with multiple inflation values
#'
#' N0sDir <- system.file("extdata", "N0-1dot3-6", package = "GeneDiscoveR")
#' overallsDir <- system.file("extdata", "Comparatives-1dot3-6", package = "GeneDiscoveR")
#' dataFile <- system.file("extdata", "annotatedCDSs.tsv", package = "GeneDiscoveR")
#' minInflation <- 1.3
#' maxInflation <- 6
#' stepInflation <- 0.1
#'
#' GeneDiscoveRobject <- GeneDiscoveR(overallsDir = overallsDir, N0sDir = N0sDir, dataFile = dataFile, minInflation = minInflation, maxInflation = maxInflation, stepInflation = stepInflation)
#'
#' # Create a GeneDiscoveR for an execution with only one inflation value
#'
#' uniqueInflation <- 1.8
#'
#' GeneDiscoveRobject <- GeneDiscoveR(N0sDir = N0sDir, dataFile = dataFile, uniqueInflation = uniqueInflation)
#'
#' # If you want to indicate the annotation file, you can use the following code:
#'
#' annotationFile <- system.file("extdata", "annotatedCDSs.tsv")
#'
#' GeneDiscoveRobject <- GeneDiscoveR(overallsDir = overallsDir, N0sDir = N0sDir, annotationFile = annotationFile, dataFile = dataFile, minInflation = minInflation, maxInflation = maxInflation, stepInflation = stepInflation)
#' @import dplyr readr
#' @importFrom GeneDiscoveR set_run_active
GeneDiscoveR <- function(overallsDir = NULL, N0sDir = NULL, dataFile = NULL, annotationFile = NULL, uniqueInflation = NULL, minInflation = NULL, maxInflation = NULL, stepInflation = NULL, pairSpeciesDir = NULL, principalSpecie = NULL, principalSpeciePrefix = NULL, orthologsTool = "OrthoFinder") {
    GeneDiscoveRobject <- list()
    if ("OrthoFinder" == orthologsTool) {
        if (is.null(N0sDir) || is.null(dataFile)) {
            stop("Error: 'overallsDir' and 'N0sDir' cannot be NULL.")
        }
        GeneDiscoveRobject$OrthologsTool <- "OrthoFinder"
        GeneDiscoveRobject$overallsDir <- overallsDir
        GeneDiscoveRobject$N0sDir <- N0sDir
        GeneDiscoveRobject$genomesData <- suppressMessages(read_tsv(file = dataFile) %>% mutate(across(everything(), as.character)))
        GeneDiscoveRobject$Nrun <- length(list.files(GeneDiscoveRobject$N0sDir, pattern = "^N0"))
        GeneDiscoveRobject$AnnotationFile <- annotationFile
        GeneDiscoveRobject$FilteredGenes <- NULL
        GeneDiscoveRobject$Phenotypes <- NULL
        GeneDiscoveRobject$RunActive <- NULL
        GeneDiscoveRobject$Identification <- NULL
        if (is.null(uniqueInflation)) {
            if (is.null(minInflation) || is.null(maxInflation) || is.null(stepInflation)) {
                stop("Error: 'minInflation', 'maxInflation', and 'stepInflation' cannot be NULL.")
            }
            GeneDiscoveRobject$Inflation <- c(minInflation, maxInflation, stepInflation)
        } else {
            GeneDiscoveRobject$Inflation <- uniqueInflation
            GeneDiscoveRobject <- set_run_active(GeneDiscoveRobject, InflationValue = uniqueInflation)
        }
    } else if ("InParanoiDB" == orthologsTool) {
        if (is.null(pairSpeciesDir) || is.null(principalSpecie) || is.null(principalSpeciePrefix)) {
            stop("Error: 'pairSpeciesDir', 'principalSpecie', and 'principalSpeciePrefix' cannot be NULL.")
        }
        GeneDiscoveRobject$OrthologsTool <- "InParanoiDB"
        GeneDiscoveRobject$pairSpeciesDir <- pairSpeciesDir
        GeneDiscoveRobject$principalSpecie <- principalSpecie
        GeneDiscoveRobject$principalSpeciePrefix <- principalSpeciePrefix
        GeneDiscoveRobject$genomesData <- suppressMessages(read_tsv(file = dataFile) %>% mutate(across(everything(), as.character)))
        GeneDiscoveRobject$Nrun <- length(list.files(GeneDiscoveRobject$pairSpeciesDir, pattern = "^inParanoid_"))
    }
    class(GeneDiscoveRobject) <- "GeneDiscoveR"

    return(GeneDiscoveRobject)
}

#' GeneDiscoveRIdentification Function
#'
#' This function creates a GeneDiscoveRIdentification object with the specified parameters.
#'
#' @param name A character string specifying the name of the object.
#' @param statistic A character string specifying the statistic to be used.
#' @param formula A formula specifying the model formula.
#' @param Phenotypes A character vector specifying the phenotypes.
#' @param columns A character vector specifying the columns to be used.
#'
#' @return A GeneDiscoveRIdentification object with the specified parameters.
#' @keywords internal
#' @examples
#' \dontrun{
#' GeneDiscoveRIdentification(name = "MyIdentification", statistic = "Fisher", formula = as.formula("A ~ B"), Phenotypes = c("A", "B"), columns = c("spIDs", "Phenotype1"))
#' }
#' @export
.GeneDiscoveRIdentification <- function(name = NULL, statistic = NULL, formula = NULL, Phenotypes = NULL, columns = NULL) {
    if (is.null(name) || is.null(statistic) || is.null(formula) || is.null(Phenotypes) || is.null(columns)) {
        stop("Error: All parameters must be provided.")
    }
    GeneDiscoveRIdentificationObject <- list()
    GeneDiscoveRIdentificationObject$name <- name
    GeneDiscoveRIdentificationObject$statistic <- statistic
    GeneDiscoveRIdentificationObject$formula <- formula
    GeneDiscoveRIdentificationObject$Phenotypes <- Phenotypes
    GeneDiscoveRIdentificationObject$columns <- columns
    class(GeneDiscoveRIdentificationObject) <- "GeneDiscoveRIdentification"
    return(GeneDiscoveRIdentificationObject)
}

#' Create a GeneDiscoveRFilteredGenes object
#'
#' This function creates a GeneDiscoveRFilteredGenes object, which is used to store filtered genes
#' based on specified criteria such as table, name, p-value, odds ratio, and sign.
#'
#' @param table The table containing the filtered genes.
#' @param name The name of the filtered genes.
#' @param pvalue The p-value threshold for filtering genes (always <= pvalue filter).
#' @param oddsRatio The odds ratio threshold for filtering genes.
#' @param sign The sign of the odds ratio for filtering genes (>, >=, <, <=).
#'
#' @return A GeneDiscoveRFilteredGenes object containing the filtered genes.
#' @keywords internal
#' @examples
#' # Create a GeneDiscoveRFilteredGenes object
#' \dontrun{
#' filtered_genes <- .GeneDiscoveRFilteredGenes(table = my_table, name = "Filtered Genes", pvalue = 0.05, oddsRatio = 2, sign = ">")
#' }
#' @export
.GeneDiscoveRFilteredGenes <- function(table = NULL, name = NULL, pvalue = NULL, oddsRatio = NULL, sign = NULL) {
    if (is.null(table) || is.null(name) || is.null(pvalue) || is.null(oddsRatio) || is.null(sign)) {
        stop("Error: All parameters must be provided.")
    }
    GeneDiscoveRFilteredGenesObject <- list()
    GeneDiscoveRFilteredGenesObject$table <- table
    GeneDiscoveRFilteredGenesObject$name <- name
    GeneDiscoveRFilteredGenesObject$pvalue <- pvalue
    GeneDiscoveRFilteredGenesObject$oddsRatio <- oddsRatio
    GeneDiscoveRFilteredGenesObject$sign <- sign
    class(GeneDiscoveRFilteredGenesObject) <- "GeneDiscoveRFilteredGenes"
    return(GeneDiscoveRFilteredGenesObject)
}

#' Selects the index of a filtered gene based on specified criteria.
#'
#' This function takes a GeneDiscoveRobject and criteria for filtering genes, such as name, p-value, odds ratio, and sign.
#' It searches through the filtered genes in the GeneDiscoveRobject and returns the index of the first gene that matches all the specified criteria.
#' If no matching gene is found, it returns NULL.
#'
#' @param GeneDiscoveRobject A GeneDiscoveR object containing filtered genes.
#' @param name The name of the identification to match.
#' @param pvalue The p-value of the filtering performed to match.
#' @param oddsRatio The odds ratio of the filtering performed to match.
#' @param sign The sign of the filtering performed to match.
#'
#' @return The index of the matching filtering performed, or NULL if no matching with any filtering performed.
#' @export
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
#  Select the index of a filtered gene
#' select_filtered_gene_index(GeneDiscoveRobject, pvalue = 0.05, oddsRatio = 1, sign = ">=", name = "PerOBtype")
#' # Output: 1
#' select_filtered_gene_index(GeneDiscoveRobject, pvalue = 0.05, oddsRatio = 1, sign = "<=", name = "PerOBtype")
#' # Output: 2
#' select_filtered_gene_index(GeneDiscoveRobject, pvalue = 0.05, oddsRatio = 1, sign = "<=", name = "OB")
#' # Output: NULL
select_filtered_gene_index <- function(GeneDiscoveRobject = NULL, name = NULL, pvalue = NULL, oddsRatio = NULL, sign = NULL) {
    if (is.null(GeneDiscoveRobject) || is.null(name) || is.null(pvalue) || is.null(oddsRatio) || is.null(sign)) {
        stop("Error: 'GeneDiscoveRobject', 'name', 'pvalue', 'oddsRatio', and 'sign' cannot be NULL.")
    }

    if (is.null(GeneDiscoveRobject$FilteredGenes)) {
        return(NULL)
    }

    for (i in seq_along(GeneDiscoveRobject$FilteredGenes)) {
        filteredGene <- GeneDiscoveRobject$FilteredGenes[[i]]
        if (filteredGene$name == name && filteredGene$pvalue == pvalue && filteredGene$oddsRatio == oddsRatio && filteredGene$sign == sign) {
            return(i)
        }
    }

    return(NULL)
}

#' Set the active run for GeneDiscoveR
#'
#' This function sets the active run for the GeneDiscoveR object by updating the relevant fields in the object.
#'
#' @param GeneDiscoveRobject The GeneDiscoveR object for which the active run needs to be set.
#' @param InflationValue The inflation value to be set for the active run, only for OrthoFinder run.
#' @param cores The number of cores to be used for parallel processing (default is 1), only for InParanoiDB run.
#'
#' @return The updated GeneDiscoveR object with the active run set.
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
#' @export
set_run_active <- function(GeneDiscoveRobject = NULL, InflationValue = 1.8, cores = 1) {
    if (is.null(GeneDiscoveRobject)) {
        stop("Error: 'GeneDiscoveRobject' cannot be NULL.")
    }

    if (GeneDiscoveRobject$OrthologsTool == "OrthoFinder") {
        if (length((GeneDiscoveRobject$Inflation)) > 1) {
            if (InflationValue < GeneDiscoveRobject$Inflation[1] || InflationValue > GeneDiscoveRobject$Inflation[2]) {
                stop("Error: 'InflationValue' is out of range.")
            }
        }
        cat("-----------From OrthoFinder-----------\n")
        GeneDiscoveRobject$RunActive$InflationActive <- InflationValue

        if (length(GeneDiscoveRobject$Inflation) > 1) {
            GeneDiscoveRobject$RunActive$IndexActive <- .index_N0_inflation(GeneDiscoveRobject, InflationValue)
        } else {
            GeneDiscoveRobject$RunActive$IndexActive <- 1
        }
        GeneDiscoveRobject$RunActive$N0fileActive <- list.files(path = GeneDiscoveRobject$N0sDir, pattern = "^N0", full.names = TRUE)[GeneDiscoveRobject$RunActive$IndexActive]
        GeneDiscoveRobject$RunActive$N0Active <- suppressMessages(read_tsv(file = GeneDiscoveRobject$RunActive$N0fileActive))
        cat("The process has been completed successfully")
    } else if (GeneDiscoveRobject$OrthologsTool == "InParanoiDB") {
        # InparanoidDB
        cat("-----------From InparanoidDB-----------\n")
        cat("-----------Step 1 of 5 - Merge inParanoid files-----------\n")
        dfMerge <- .merge_inParanoid_files(GeneDiscoveRobject$pairSpeciesDir)
        cat("-----------Step 2 of 5 - Merge inParanoid files-----------\n")
        results <- .check_OrtoA_in_OrtoB(dfMerge, cores = cores)
        cat("-----------Step 3 of 5 - Merge inParanoid files-----------\n")
        dfMerge <- .merge_groups(dfMerge, results)
        cat("-----------Step 4 of 5 - Merge inParanoid files-----------\n")
        dfMerge$dfMergemodifyunprocess <- .merge_groups_parallel(dfMerge$dfMergemodifyunprocess, cores = cores)
        dfMerge <- rbind(dfMerge$dfMergemodifyunprocess, dfMerge$dfMergemodifyprocess)
        dfMerge <- dfMerge %>% rename(OrtoA = "dfMerge$OrtoA[i]")
        cat("-----------Step 5 of 5 - Merge inParanoid files-----------\n")
        dfMerge <- .update_dfMerge(dfMerge, GeneDiscoveRobject$principalSpeciePrefix, GeneDiscoveRobject$genomesData, colnames(GeneDiscoveRobject$genomesData)[1])
        cat("The process has been completed successfully\n\n")
        GeneDiscoveRobject$RunActive$N0Active <- dfMerge
        GeneDiscoveRobject$RunActive$N0Active <- GeneDiscoveRobject$RunActive$N0Active %>% mutate(OG = 1:nrow(GeneDiscoveRobject$RunActive$N0Active), HOG = 1:nrow(GeneDiscoveRobject$RunActive$N0Active))
        GeneDiscoveRobject$RunActive$N0Active
    }
    return(GeneDiscoveRobject)
}

#' Clean phenotypes
#'
#' This function removes the phenotypes from a GeneDiscoveR object.
#'
#' @param GeneDiscoveRobject A GeneDiscoveR object.
#' @return The GeneDiscoveR object with the phenotypes removed.
#' @export
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
#' # Phenotypes
#' names(GeneDiscoveRobject$Phenotypes)
#' # Output: [1] "one_in_specialized_cell" "many_in_all_cells"
#'
#' # Clean phenotypes
#' GeneDiscoveRobject <- clean_phenotypes(GeneDiscoveRobject)
#' names(GeneDiscoveRobject$Phenotypes)
#' # Output: NULL
clean_phenotypes <- function(GeneDiscoveRobject = NULL) {
    if (is.null(GeneDiscoveRobject)) {
        stop("Error: 'GeneDiscoveRobject' cannot be NULL.")
    }

    GeneDiscoveRobject$Phenotypes <- NULL

    return(GeneDiscoveRobject)
}

#' Calculate the index of a specific inflation value in a GeneDiscoveR object.
#'
#' This function calculates the index of a specific inflation value in a GeneDiscoveR object.
#'
#' @param GeneDiscoveRobject A GeneDiscoveR object.
#' @param InflationValue The inflation value for which the index needs to be calculated.
#'
#' @return The index of the specified inflation value in the GeneDiscoveR object.
#' @keywords internal
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
#' # Calculate the index of inflation value 1.8
#' index <- .index_N0_inflation(GeneDiscoveRobject, 1.8)
#' index
#' # Output: 6
#' @export
.index_N0_inflation <- function(GeneDiscoveRobject = NULL, InflationValue = 1.8) {
    if (is.null(GeneDiscoveRobject)) {
        stop("Error: 'GeneDiscoveRobject' cannot be NULL.")
    }

    index <- which(seq(GeneDiscoveRobject$Inflation[1], GeneDiscoveRobject$Inflation[2], GeneDiscoveRobject$Inflation[3]) == InflationValue)

    return(index)
}

#' Calculate the sum per phenotype
#'
#' This function calculates the sum of non-missing values per phenotype in a data frame.
#'
#' @param df The input data frame.
#' @param predictor The name of the predictor variable.
#' @param response The name of the response variable.
#' @param nameColumn1 The name of the column to store the sum for the predictor variable.
#' @param nameColumn2 The name of the column to store the sum for the response variable.
#'
#' @return The modified data frame with additional columns storing the sums per phenotype.
#' @keywords internal
#' @examples
#' \dontrun{
#' df <- data.frame(predictor = c(1, 2, NA, 4), response = c(NA, 2, 3, 4))
#' .sum_per_phenotype(df, "predictor", "response", "sum_predictor", "sum_response")
#' }
#' @import dplyr
#' @export
.sum_per_phenotype <- function(df, predictor, response, nameColumn1, nameColumn2) {
    if (!(!!nameColumn1 %in% colnames(df)) && !(!!nameColumn2 %in% colnames(df))) {
        df <- df %>%
            rowwise() %>%
            mutate(
                nameColumn1 = sum(!is.na(c_across(all_of(predictor))) & c_across(all_of(predictor)) != ""),
                nameColumn2 = sum(!is.na(c_across(all_of(response))) & c_across(all_of(response)) != "")
            )
    } else if (!!nameColumn1 %in% colnames(df)) {
        df <- df %>%
            rowwise() %>%
            mutate(
                nameColumn2 = sum(!is.na(c_across(all_of(response))) & c_across(all_of(response)) != "")
            ) %>%
            rename(nameColumn1 := !!nameColumn1)
    } else if (!!nameColumn2 %in% colnames(df)) {
        df <- df %>%
            rowwise() %>%
            mutate(
                nameColumn1 = sum(!is.na(c_across(all_of(predictor))) & c_across(all_of(predictor)) != "")
            ) %>%
            rename(nameColumn2 := !!nameColumn2)
    } else {
        df <- df %>% rename(nameColumn1 := !!nameColumn1, nameColumn2 := !!nameColumn2)
    }
    return(df)
}

#' Calculate Fisher's exact test for each row in a data frame
#'
#' This function calculates Fisher's exact test for each row in a data frame. It takes a data frame, a predictor variable, a response variable, and three column names as input. It performs Fisher's exact test for each row and adds the test result as a new column to the data frame.
#'
#' @param df A data frame containing the predictor and response variables
#' @param predictor A vector representing the predictor variable
#' @param response A vector representing the response variable
#' @param nameColumn1 The name of the column in the data frame representing the first variable in the contingency table
#' @param nameColumn2 The name of the column in the data frame representing the second variable in the contingency table
#' @param nameColumn3 The name of the column to be added to the data frame to store the Fisher's exact test result
#'
#' @return The input data frame with an additional column containing the Fisher's exact test result for each row
#' @keywords internal
#' @examples
#' \dontrun{
#' df <- data.frame(predictor = c(1, 0, 1), response = c(1, 1, 0))
#' df <- .fisher_per_row(df, df$predictor, df$response, "nameColumn1", "nameColumn2", "nameColumn3")
#' print(df)
#' }
#' @import dplyr
#' @export
.fisher_per_row <- function(df, predictor, response, nameColumn1, nameColumn2, nameColumn3) {
    df <- df %>%
        rowwise() %>%
        mutate(
            !!nameColumn3 := list(fisher.test(matrix(
                c(
                    nameColumn1,
                    nameColumn2,
                    length(predictor) - nameColumn1,
                    length(response) - nameColumn2
                ),
                nrow = 2, byrow = TRUE
            )))
        )
    return(df)
}

#' Get the names of the identification objects in a GeneDiscoveR object
#'
#' This function takes a GeneDiscoveR object as input and returns a vector
#' containing the names of the identification objects within the GeneDiscoveR
#' object.
#'
#' @param GeneDiscoveRobject A GeneDiscoveR object.
#'
#' @return A character vector containing the names of the identification objects.
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
#' GeneDiscoveRobject <- gene_identification_by_phenotype(GeneDiscoveRobject = GeneDiscoveRobject, formula = as.formula("one_in_specialized_cell ~ many_in_all_cells"), statistic = "Fisher", name = "OneInSpeVSmanyInAll", cores = 8)
#' GeneDiscoveRobject <- gene_identification_by_phenotype(GeneDiscoveRobject = GeneDiscoveRobject, formula = as.formula("many_in_all_cells ~ one_in_specialized_cell"), statistic = "Fisher", name = "manyInAllVSOneInSpe", cores = 8)
#'
#' # Get the names of the identification objects
#' get_names_identification(GeneDiscoveRobject)
#' # Output: [1] "OneInSpeVSmanyInAll" "manyInAllVSOneInSpe"
#'
#' @export
get_names_identification <- function(GeneDiscoveRobject = NULL) {
    if (is.null(GeneDiscoveRobject)) {
        stop("Error: 'GeneDiscoveRobject' cannot be NULL.")
    }
    if (is.null(GeneDiscoveRobject$Identification)) {
        stop("Error: 'GeneDiscoveRobject$Identification' cannot be NULL.")
    }
    result <- c()
    for (i in seq_along(GeneDiscoveRobject$Identification)) {
        result <- c(result, GeneDiscoveRobject$Identification[[i]]$name)
    }
    return(result)
}

#' Get complete table from GeneDiscoveR object
#'
#' This function retrieves the complete table from a GeneDiscoveR object.
#' If you perform a gene identification, the complete table will be the table with the results of the gene identification.
#'
#' @param GeneDiscoveRobject A GeneDiscoveR object.
#'
#' @return The complete table from the GeneDiscoveR object.
#'
#' @examples
#  Create a GeneDiscoveR object
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
#' # Get complete table
#' table <- get_complete_table(GeneDiscoveRobject)
#' # Output: A data frame with the complete table from the GeneDiscoveR object.
#' @export
get_complete_table <- function(GeneDiscoveRobject = NULL) {
    if (is.null(GeneDiscoveRobject)) {
        stop("Error: 'GeneDiscoveRobject' cannot be NULL.")
    }
    if (is.null(GeneDiscoveRobject$RunActive$N0Active)) {
        stop("Error: 'GeneDiscoveRobject$RunActive' cannot be NULL.")
    }
    excludecols <- c("original_index", "contains-gene", "log-odds-ratio")
    excludecols <- c(excludecols, grep("^fisherResult", names(GeneDiscoveRobject$RunActive$N0Active), value = TRUE))
    result <- GeneDiscoveRobject$RunActive$N0Active[, -which(names(GeneDiscoveRobject$RunActive$N0Active) %in% excludecols)]
    return(result)
}
