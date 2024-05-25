#' Run the GeneDiscover web application
#'
#' This function launches the GeneDiscover web application.
#'
#' @return None
#' @examples
#' \dontrun{
#' run_genediscover_web_app()
#' }
#' @export
run_genediscover_web_app <- function() {
    if (!requireNamespace("shiny", quietly = TRUE)) {
        install.packages("shiny")
        if (!requireNamespace("shiny", quietly = TRUE)) {
            stop("Failed to install the 'shiny' package.")
        }
    }

    if (!requireNamespace("plotly", quietly = TRUE)) {
        install.packages("plotly")
        if (!requireNamespace("plotly", quietly = TRUE)) {
            stop("Failed to install the 'plotly' package.")
        }
    }

    if (!requireNamespace("shinydashboard", quietly = TRUE)) {
        install.packages("shinydashboard")
        if (!requireNamespace("shinydashboard", quietly = TRUE)) {
            stop("Failed to install the 'shinydashboard' package.")
        }
    }

    if (!requireNamespace("DT", quietly = TRUE)) {
        install.packages("DT")
        if (!requireNamespace("DT", quietly = TRUE)) {
            stop("Failed to install the 'DT' package.")
        }
    }
    launch_genediscover_web_app()
}

#' Get identification by name
#'
#' This function retrieves the identification information from a GeneDiscoveR object
#' based on the specified name.
#'
#' @param GeneDiscoveRobject A GeneDiscoveR object.
#' @param name The name of the identification to retrieve.
#'
#' @return The identification information matching the specified name.
#' @keywords internal
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
#' # Get identification by name
#' .get_identification(GeneDiscoveRobject, "PerOBtype")
#' @export
.get_identification <- function(GeneDiscoveRobject = NULL, name = NULL) {
    if (is.null(GeneDiscoveRobject) || is.null(name)) {
        stop("Error: 'GeneDiscoveRobject' and 'name' cannot be NULL.")
    }
    if (is.null(GeneDiscoveRobject$Identification)) {
        stop("Error: 'GeneDiscoveRobject$Identification' cannot be NULL.")
    }
    for (i in seq_along(GeneDiscoveRobject$Identification)) {
        if (GeneDiscoveRobject$Identification[[i]]$name == name) {
            return(GeneDiscoveRobject$Identification[[i]])
        }
    }
}

#' Get Filtered Genes Table
#'
#' This function retrieves a filtered genes table based on specified criteria.
#'
#' @param GeneDiscoveRobject A GeneDiscoveR object containing the filtered genes.
#' @param name The name of the gene.
#' @param pvalue The p-value of the gene.
#' @param oddsRatio The odds ratio of the gene.
#' @param sign The sign of the gene.
#'
#' @return A filtered genes table that matches the specified criteria.
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
#'
#' # Get filtered genes table
#' get_filtered_genes_table(GeneDiscoveRobject, name = "PerOBtype", pvalue = 0.05, oddsRatio = 1, sign = ">=")
#' @export
get_filtered_genes_table <- function(GeneDiscoveRobject = NULL, name = NULL, pvalue = NULL, oddsRatio = NULL, sign = NULL) {
    if (is.null(GeneDiscoveRobject) || is.null(name) || is.null(pvalue) || is.null(oddsRatio) || is.null(sign)) {
        stop("Error: 'GeneDiscoveRobject', 'name', 'pvalue', 'oddsRatio', and 'sign' cannot be NULL.")
    }
    if (is.null(GeneDiscoveRobject$FilteredGenes)) {
        stop("Error: 'GeneDiscoveRobject$FilteredGenes' cannot be NULL.")
    }
    result <- c()
    for (i in seq_along(GeneDiscoveRobject$FilteredGenes)) {
        filteredGene <- GeneDiscoveRobject$FilteredGenes[[i]]
        if (filteredGene$name == name && filteredGene$pvalue == pvalue && filteredGene$oddsRatio == oddsRatio && filteredGene$sign == sign) {
            result <- filteredGene$table[, !startsWith(names(filteredGene$table), "fisherResult")]
        }
    }
    return(result)
}

#' Plot GeneDiscoveR Volcano
#'
#' This function plots a volcano plot for GeneDiscoveR analysis results.
#' With the \code{formula("first_phenotype ~ second_phenotype")},
#' The volcano plot for odds ratio > 1 represent the orthogroups enriched with genes from species with the first phenotype.
#' In contrast, the volcano plot for odds ratio < 1 represent the orthogroups enriched with genes from species with the second phenotype.
#'
#' @param GeneDiscoveRobject An object of class 'GeneDiscoveR' containing the analysis results.
#' @param name The name of the analysis to plot.
#'
#' @return A ggplot object representing the volcano plot.
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
#' # Volcano plot
#' plot <- plot_genediscover_volcano(GeneDiscoveRobject, name = "PerOBtype")
#' # Output: A ggplot object representing the volcano plot.
#' @import ggplot2 dplyr ggsci
#' @export
plot_genediscover_volcano <- function(GeneDiscoveRobject = NULL, name = NULL) {
    GeneDiscoveRidentification <- .get_identification(GeneDiscoveRobject = GeneDiscoveRobject, name = name)
    if (!inherits(GeneDiscoveRidentification, "GeneDiscoveRIdentification")) {
        stop("Error: 'GeneDiscoveRidentification' must be of class 'GeneDiscoveRIdentification'.")
    }
    table <- GeneDiscoveRobject$RunActive$N0Active %>%
        mutate(
            logoddRatioFisher = case_when(
                !!sym(GeneDiscoveRidentification$columns[5]) == 0 ~ -5,
                !!sym(GeneDiscoveRidentification$columns[5]) == Inf ~ 5,
                is.finite(!!sym(GeneDiscoveRidentification$columns[5])) ~ log(!!sym(GeneDiscoveRidentification$columns[5]))
            )
        )
    .set_ggplot2_theme()
    g1 <- ggplot() +
        geom_point(aes(x = -log(table[[GeneDiscoveRidentification$columns[4]]]), y = table$logoddRatioFisher), color = "#696D7D") +
        coord_flip() +
        geom_hline(yintercept = 0, linetype = "dotted", col = "red") +
        geom_vline(xintercept = 2.9957, linetype = "dotted", col = "red") +
        annotate("text", x = 2.9957, y = 0, label = "p-value <= 0.05", vjust = 1.5, color = "black") +
        annotate("text", x = 8, y = 1, label = "OddRatio >= 1", color = "black") +
        ylab("log(Odds Ratio)") +
        xlab("-log(p-value)")
    return(g1)
}

#' Get overall statistics from GeneDiscoveR object
#'
#' This function retrieves the overall statistics from a GeneDiscoveR object.
#'
#' @param GeneDiscoveRobject A GeneDiscoveR object.
#'
#' @return A data frame containing the overall metrics from the GeneDiscoveR object.
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
#' # Calculate overall statistics
#'
#' GeneDiscoveRobject <- calculate_overall_statistics(GeneDiscoveRobject, cores = 8)
#'
#' # Get overall statistics
#' overallStatistics <- get_overall_statistics(GeneDiscoveRobject)
#' # Output: A data frame containing the overall metrics from the GeneDiscoveR object.
#'
#' @seealso
#' \code{\link{GeneDiscoveR}}
#' \code{\link{calculate_overall_statistics}}
#' \code{\link{plot_allSpeciesOGs_sOGs_per_inflation}}
get_overall_statistics <- function(GeneDiscoveRobject = NULL) {
    if (is.null(GeneDiscoveRobject)) {
        stop("Error: 'GeneDiscoveRobject' cannot be NULL.")
    }
    if (is.null(GeneDiscoveRobject$runsData$overallMetrics)) {
        stop("Error: 'GeneDiscoveRobject$runsData$overallMetrics' cannot be NULL.")
    }
    return(GeneDiscoveRobject$runsData$overallMetrics)
}

#' Get the median statistics from a GeneDiscoveR object
#'
#' This function retrieves the median statistics from a GeneDiscoveR object.
#'
#' @param GeneDiscoveRobject A GeneDiscoveR object.
#' @return A numeric vector containing the median statistics.
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
#' # Calculate median statistics
#' GeneDiscoveRobject <- calculate_median_statistics(GeneDiscoveRobject, cores = 8)
#'
#' # Get median statistics
#' medianStatistics <- get_medians_statistics(GeneDiscoveRobject)
#' # Output: A data frame containing the median statistics from the GeneDiscoveR object.
#' @seealso
#' \code{\link{GeneDiscoveR}}
#' \code{\link{calculate_median_statistics}}
#' \code{\link{plot_OGs_HOGs_per_inflation}}
#' @export
get_medians_statistics <- function(GeneDiscoveRobject = NULL) {
    if (is.null(GeneDiscoveRobject)) {
        stop("Error: 'GeneDiscoveRobject' cannot be NULL.")
    }
    if (is.null(GeneDiscoveRobject$runsData$medians)) {
        stop("Error: 'GeneDiscoveRobject$runsData$medians' cannot be NULL.")
    }
    return(GeneDiscoveRobject$runsData$medians)
}

#' Set ggplot2 theme
#'
#' This function sets a custom theme for ggplot2 plots.
#'
#' @details The function sets various theme elements for ggplot2 plots, including the legend background,
#' plot title and subtitle, panel grid and background, axis text and title, legend key, and panel border.
#' The function also sets the theme using `theme_set()` function from ggplot2 package.
#'
#' @import ggplot2
#'
#' @return This function does not return any value.
#' @keywords internal
#' @examples
#' .set_ggplot2_theme()
#' @export
.set_ggplot2_theme <- function() {
    theme <- theme(
        legend.background = element_blank(),
        plot.title = element_text(
            hjust = 0.5,
            vjust = 2,
            size = 16
        ),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = -90),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(vjust = 1, size = 15, angle = 90),
        legend.key = element_blank(),
        panel.border = element_rect(
            colour = "#77767b",
            fill = NA,
            linewidth = 0.5
        )
    )
    theme_set(theme)
}
