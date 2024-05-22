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

#' Create a GeneDiscoveR object
#'
#' This function creates a GeneDiscoveR object with the specified directories for overalls and N0s.
#'
#' @param overallsDir The directory path for the overalls.
#' @param N0sDir The directory path for the N0s.
#'
#' @return A GeneDiscoveR object with the overalls directory, N0s directory, and the number of runs.
#' @export
#'
#' @examples
#' # Create a GeneDiscoveR object with the specified directories
#' pheno_obj <- create_GeneDiscoveR_object(overallsDir = "/path/to/overalls", N0sDir = "/path/to/N0s")
#'
#' # Access the overalls directory
#' pheno_obj$overallsDir
#'
#' # Access the N0s directory
#' pheno_obj$N0sDir
#'
#' # Access the number of runs
#' pheno_obj$Nrun
GeneDiscoveR <- function(overallsDir = NULL, N0sDir = NULL, dataFile = NULL, annotationFile = NULL, uniqueInflation = NULL, minInflation = NULL, maxInflation = NULL, stepInflation = NULL) {
    if (is.null(N0sDir) || is.null(dataFile)) {
        stop("Error: 'overallsDir' and 'N0sDir' cannot be NULL.")
    }

    GeneDiscoveRobject <- list()
    GeneDiscoveRobject$overallsDir <- overallsDir
    GeneDiscoveRobject$N0sDir <- N0sDir
    GeneDiscoveRobject$genomesData <- read_tsv(file = dataFile) %>% mutate(across(everything(), as.character))
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

    class(GeneDiscoveRobject) <- "GeneDiscoveR"

    return(GeneDiscoveRobject)
}

GeneDiscoveRIdentification <- function(name = NULL, statistic = NULL, formula = NULL, Phenotypes = NULL, columns = NULL) {
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

GeneDiscoveRFilteredGenes <- function(table = NULL, name = NULL, pvalue = NULL, oddsRatio = NULL, sign = NULL) {
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

set_run_active <- function(GeneDiscoveRobject = NULL, InflationValue = 1.8) {
    if (is.null(GeneDiscoveRobject)) {
        stop("Error: 'GeneDiscoveRobject' cannot be NULL.")
    }

    if (length((GeneDiscoveRobject$Inflation)) > 1) {
        if (InflationValue < GeneDiscoveRobject$InflationLimits[1] || InflationValue > GeneDiscoveRobject$InflationLimits[2]) {
            stop("Error: 'InflationValue' is out of range.")
        }
    }

    GeneDiscoveRobject$RunActive$InflationActive <- InflationValue

    if (length(GeneDiscoveRobject$Inflation) > 1) {
        GeneDiscoveRobject$RunActive$IndexActive <- index_N0_inflation(GeneDiscoveRobject, InflationValue)
    } else {
        GeneDiscoveRobject$RunActive$IndexActive <- 1
    }
    GeneDiscoveRobject$RunActive$N0fileActive <- list.files(path = GeneDiscoveRobject$N0sDir, pattern = "^N0", full.names = TRUE)[GeneDiscoveRobject$RunActive$IndexActive]
    GeneDiscoveRobject$RunActive$N0Active <- read_tsv(file = GeneDiscoveRobject$RunActive$N0fileActive)

    return(GeneDiscoveRobject)
}

clean_phenotypes <- function(GeneDiscoveRobject = NULL) {
    if (is.null(GeneDiscoveRobject)) {
        stop("Error: 'GeneDiscoveRobject' cannot be NULL.")
    }

    GeneDiscoveRobject$Phenotypes <- NULL

    return(GeneDiscoveRobject)
}

index_N0_inflation <- function(GeneDiscoveRobject = NULL, InflationValue = 1.8) {
    if (is.null(GeneDiscoveRobject)) {
        stop("Error: 'GeneDiscoveRobject' cannot be NULL.")
    }

    index <- which(seq(GeneDiscoveRobject$InflationLimits[1], GeneDiscoveRobject$InflationLimits[2], GeneDiscoveRobject$InflationLimits[3]) == InflationValue)

    return(index)
}

process_rows <- function(df, predictor, response, nameColumn1, nameColumn2) {
    df <- df %>%
        rowwise() %>%
        mutate(
            nameColumn1 = sum(!is.na(c_across(all_of(predictor)))),
            nameColumn2 = sum(!is.na(c_across(all_of(response))))
        )
    return(df)
}

process_rows2 <- function(df, predictor, response, nameColumn1, nameColumn2, nameColumn3) {
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
