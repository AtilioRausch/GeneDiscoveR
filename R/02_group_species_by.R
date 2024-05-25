#' Select species by phenotype
#'
#' This function selects species from a GeneDiscoveR object based on a specified phenotype.
#'
#' @param GeneDiscoveRobject A GeneDiscoveR object containing genome data.
#' @param columnPhenotype The name of the column in the genome data that contains the phenotype information.
#' @param columnID The name of the column in the genome data that contains the species ID.
#' @param type The specific phenotype to filter by.
#'
#' @return The updated GeneDiscoveR object with the selected species.
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
#' GeneDiscoveRobject <- select_species_by_phenotype(GeneDiscoveRobject = GeneDiscoveRobject, columnPhenotype = "Oil-body-type", columnID = "OrthofinderID", type = "noneOB")
#'
#' # Show Phenotypes
#' print(GeneDiscoveRobject$Phenotypes)
#' # Output: $one_in_specialized_cell ... $many_in_all_cells ... $noneOB ...
#'
#' # Clean phenotypes
#' GeneDiscoveRobject <- clean_phenotypes(GeneDiscoveRobject)
#' print(GeneDiscoveRobject$Phenotypes)
#' # Output: NULL
#' @export
#' @import dplyr
select_species_by_phenotype <- function(GeneDiscoveRobject = NULL, columnPhenotype = "Oil-body-type", columnID = "OrthofinderID", type = NULL) {
    if (is.null(GeneDiscoveRobject) || is.null(GeneDiscoveRobject$genomesData)) {
        stop("Please provide the GeneDiscoveRobject parameter and make sure GeneDiscoveRobject$genomesData is not null.")
    }
    result <- GeneDiscoveRobject$genomesData %>%
        filter(!!sym(columnPhenotype) == type) %>%
        select(!!sym(columnID))
    result <- .table_as_character(table = result)

    GeneDiscoveRobject$Phenotypes[[type]] <- result

    return(GeneDiscoveRobject)
}
