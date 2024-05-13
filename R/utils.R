#' Convert a table to character vector
#'
#' This function takes a table as input and converts it to a character vector.
#'
#' @param table The table to be converted.
#' @return A character vector representing the table.
#' @examples
#' table <- data.frame(x = c(1, 2, 3), y = c("a", "b", "c"))
#' table_as_character(table)
#' # Output: [1] "1" "2" "3" "a" "b" "c"
#'
#' @import dplyr
#' @export
table_as_character <- function(table = NULL) {
    if (is.null(table)) {
        stop("Error: 'table' cannot be NULL.")
    }

    result <- table %>%
        pull(.) %>%
        as.character()

    return(result)
}

#' Create a PhenoR object
#'
#' This function creates a PhenoR object with the specified directories for overalls and N0s.
#'
#' @param overallsDir The directory path for the overalls.
#' @param N0sDir The directory path for the N0s.
#'
#' @return A PhenoR object with the overalls directory, N0s directory, and the number of runs.
#' @export
#'
#' @examples
#' # Create a PhenoR object with the specified directories
#' pheno_obj <- create_phenor_object(overallsDir = "/path/to/overalls", N0sDir = "/path/to/N0s")
#'
#' # Access the overalls directory
#' pheno_obj$overallsDir
#'
#' # Access the N0s directory
#' pheno_obj$N0sDir
#'
#' # Access the number of runs
#' pheno_obj$Nrun
PhenoR <- function(overallsDir = NULL, N0sDir = NULL, dataFile = NULL, annotationFile = NULL, minInflation = 1.3, maxInflation = 6, stepInflation = 0.1) {
    if (is.null(overallsDir) || is.null(N0sDir) || is.null(dataFile)) {
        stop("Error: 'overallsDir' and 'N0sDir' cannot be NULL.")
    }
    PhenoRobject <- list()
    PhenoRobject$overallsDir <- overallsDir
    PhenoRobject$N0sDir <- N0sDir
    PhenoRobject$genomesData <- read_tsv(file = dataFile)
    PhenoRobject$Nrun <- length(list.files(PhenoRobject$N0sDir))
    PhenoRobject$InflationLimits <- c(minInflation, maxInflation, stepInflation)
    PhenoRobject$AnnotationFile <- annotationFile
    PhenoRobject$FilteredGenes <- NULL
    PhenoRobject$Phenotypes <- NULL
    PhenoRobject$RunActive <- NULL
    PhenoRobject$Identification <- NULL

    class(PhenoRobject) <- "PhenoR"

    return(PhenoRobject)
}

PhenoRIdentification <- function(name = NULL, statistic = NULL, formula = NULL, Phenotypes = NULL, columns = NULL) {
    if (is.null(name) || is.null(statistic) || is.null(formula) || is.null(Phenotypes) || is.null(columns)) {
        stop("Error: All parameters must be provided.")
    }
    PhenoRIdentificationObject <- list()
    PhenoRIdentificationObject$name <- name
    PhenoRIdentificationObject$statistic <- statistic
    PhenoRIdentificationObject$formula <- formula
    PhenoRIdentificationObject$Phenotypes <- Phenotypes
    PhenoRIdentificationObject$columns <- columns
    class(PhenoRIdentificationObject) <- "PhenoRIdentification"
    return(PhenoRIdentificationObject)
}
PhenoRFilteredGenes <- function(table = NULL, name = NULL, pvalue = NULL, oddsRatio = NULL, sign = NULL) {
    if (is.null(table) || is.null(name) || is.null(pvalue) || is.null(oddsRatio) || is.null(sign)) {
        stop("Error: All parameters must be provided.")
    }
    PhenoRFilteredGenesObject <- list()
    PhenoRFilteredGenesObject$table <- table
    PhenoRFilteredGenesObject$name <- name
    PhenoRFilteredGenesObject$pvalue <- pvalue
    PhenoRFilteredGenesObject$oddsRatio <- oddsRatio
    PhenoRFilteredGenesObject$sign <- sign
    class(PhenoRFilteredGenesObject) <- "PhenoRFilteredGenes"
    return(PhenoRFilteredGenesObject)
}

select_filtered_gene_index <- function(PhenoRobject = NULL, name = NULL, pvalue = NULL, oddsRatio = NULL, sign = NULL) {
    if (is.null(PhenoRobject) || is.null(name) || is.null(pvalue) || is.null(oddsRatio) || is.null(sign)) {
        stop("Error: 'PhenoRobject', 'name', 'pvalue', 'oddsRatio', and 'sign' cannot be NULL.")
    }

    if (is.null(PhenoRobject$FilteredGenes)) {
        return(NULL)
    }

    for (i in seq_along(PhenoRobject$FilteredGenes)) {
        filteredGene <- PhenoRobject$FilteredGenes[[i]]
        if (filteredGene$name == name && filteredGene$pvalue == pvalue && filteredGene$oddsRatio == oddsRatio && filteredGene$sign == sign) {
            return(i)
        }
    }

    return(NULL)
}

set_run_active <- function(PhenoRobject = NULL, InflationValue = 1.8) {
    if (is.null(PhenoRobject)) {
        stop("Error: 'PhenoRobject' cannot be NULL.")
    }

    if (InflationValue < PhenoRobject$InflationLimits[1] || InflationValue > PhenoRobject$InflationLimits[2]) {
        stop("Error: 'InflationValue' is out of range.")
    }

    PhenoRobject$RunActive$InflationActive <- InflationValue
    PhenoRobject$RunActive$IndexActive <- index_N0_inflation(PhenoRobject, InflationValue)
    PhenoRobject$RunActive$N0fileActive <- list.files(path = PhenoRobject$N0sDir, pattern = "^N0", full.names = TRUE)[PhenoRobject$RunActive$IndexActive]
    PhenoRobject$RunActive$N0Active <- read_tsv(file = PhenoRobject$RunActive$N0fileActive)

    return(PhenoRobject)
}

clean_phenotypes <- function(PhenoRobject = NULL) {
    if (is.null(PhenoRobject)) {
        stop("Error: 'PhenoRobject' cannot be NULL.")
    }

    PhenoRobject$Phenotypes <- NULL

    return(PhenoRobject)
}

index_N0_inflation <- function(PhenoRobject = NULL, InflationValue = 1.8) {
    if (is.null(PhenoRobject)) {
        stop("Error: 'PhenoRobject' cannot be NULL.")
    }

    index <- which(seq(PhenoRobject$InflationLimits[1], PhenoRobject$InflationLimits[2], PhenoRobject$InflationLimits[3]) == InflationValue)

    return(index)
}
