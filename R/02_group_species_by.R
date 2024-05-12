#' Select species by oil body type
#'
#' This function selects species based on their oil body type from a set of genomes and transcriptomes.
#'
#' @param genomes A data frame containing information about genomes.
#' @param transcriptomes A data frame containing information about transcriptomes.
#' @param columnPhenotype The name of the column in the data frames that represents the oil body type. Default is `Oil-body-type`.
#' @param columnOrthofinderID The name of the column in the data frames that represents the Orthofinder ID. Default is `OrthofinderID`.
#' @param type The oil body type to filter by. Default is NULL, which selects all types.
#' @param origin The origin of the data to filter. Possible values are "G&T" (genomes and transcriptomes), "G" (genomes only), and "T" (transcriptomes only).
#'
#' @return A table with the selected species' Orthofinder names.
#'
#' @examples
#' genomes <- data.frame(`Oil-body-type` = c("Type A", "Type B", "Type C"), OrthofinderID = c("Species1", "Species2", "Species3"))
#' transcriptomes <- data.frame(`Oil-body-type` = c("Type A", "Type B", "Type C"), OrthofinderID = c("Species4", "Species5", "Species6"))
#' select_species_by_phenotype(genomes, transcriptomes, type = "Type A", origin = "G&T")
#'
#' @import dplyr
#' @import tidyr
#' @export
select_species_by_phenotype <- function(PhenoRobject = NULL, columnPhenotype = `Oil-body-type`, columnOrthofinderID = `OrthofinderID`, type = NULL) {
    if (is.null(PhenoRobject) || is.null(PhenoRobject$genomesData)) {
        stop("Please provide the PhenoRobject parameter and make sure PhenoRobject$genomesData is not null.")
    }
    result <- PhenoRobject$genomesData %>%
        filter({{ columnPhenotype }} == type) %>%
        select({{ columnOrthofinderID }})
    result <- table_as_character(table = result)

    PhenoRobject$Phenotypes[[type]] <- result

    return(PhenoRobject)
}
pnameColumn3 <- paste0("fisherResult", name)
            PhenoRobject$RunActive$N0Active <- PhenoRobject$RunActive$N0Active %>%
                rowwise() %>%
                mutate(
                    !!nameColumn3 := list(fisher.test(matrix(
                        c(
                            !!sym(nameColumn1),
                            !!sym(nameColumn2)
                        ),
                        nrow = 2, byrow = TRUE
                    ),
                    length(PhenoRobject$Phenotypes[[formula[[2]]]]) - !!sym(nameColumn1),
                    length(PhenoRobject$Phenotypes[[formula[[3]]]]) - !!sym(nameColumn2)
                    ))),
                    !!paste0("pvalueFisher", name) := nameColumn3[[1]],
                    !!paste0("oddsRatioFisher", name) := nameColumn3[[2]]
                )objectPhenoRobject