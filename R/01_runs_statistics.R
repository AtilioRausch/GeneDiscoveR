#' Calculate statistics and return metrics dataframe
#'
#' This function calculates obtain the overall statistics for the runs and returns the metrics dataframe.
#' The metrics include the mean, median, G50 (assigned genes), G50 (all genes), O50 (assigned genes),
#' O50 (all genes), All species OG, and sOGs per run.
#'
#' @param GeneDiscoveRobject The GeneDiscoveR object containing the overall directory and number of runs.
#' @param cores The number of cores to use for parallel processing. Default is 1.
#' @return A GeneDiscoveR object with the calculated statistics added to the runsData$overallMetrics field.
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
#' @import dplyr tidyr readr
#' @import doParallel
#' @import parallel
#' @export
calculate_overall_statistics <- function(GeneDiscoveRobject = NULL, cores = 1) {
    if (is.null(GeneDiscoveRobject) || is.null(GeneDiscoveRobject$overallsDir)) {
        stop("Please provide the GeneDiscoveRobject parameter and make sure GeneDiscoveRobject$overallsDir is not null.")
    }

    if (GeneDiscoveRobject$Nrun <= 9) {
        indexs <- paste0("0", seq(1, GeneDiscoveRobject$Nrun))
    } else {
        indexs <- paste0("0", seq(1, 9))
        indexs <- c(indexs, seq(10, GeneDiscoveRobject$Nrun))
    }

    metrics <- NULL

    # Parallelize the loop using 'foreach' and 'doParallel' packages
    if (cores > 1) {
        # Check if the number of threads is greater than the number of available cores
        if (cores > detectCores()) {
            stop("The number of cores cannot be greater than the number of available cores.")
        }

        # Create a cluster with the specified number of cores
        cl <- makeCluster(cores)
        registerDoParallel(cl)

        # Loop through the files in parallel
        metrics <- foreach(index = indexs, .combine = rbind, .packages = c("readr")) %dopar% {
            stats <- paste0(GeneDiscoveRobject$overallsDir, "/Statistics_Overall_", index, ".tsv")
            stats <- suppressMessages(readr::read_tsv(stats, col_names = FALSE, skip = 10, n_max = 8))
            cbind(
                stats[1, 2], stats[2, 2], stats[3, 2], stats[4, 2],
                stats[5, 2], stats[6, 2], stats[7, 2], stats[8, 2]
            )
        }

        # Stop the cluster
        stopCluster(cl)
    } else {
        # Loop through the files sequentially
        for (index in indexs) {
            stats <- paste0(GeneDiscoveRobject$overallsDir, "/Statistics_Overall_", index, ".tsv")
            stats <- suppressMessages(readr::read_tsv(stats, col_names = FALSE, skip = 10, n_max = 8))
            metrics <- rbind(metrics, cbind(
                stats[1, 2], stats[2, 2], stats[3, 2], stats[4, 2],
                stats[5, 2], stats[6, 2], stats[7, 2], stats[8, 2]
            ))
        }
    }

    metrics <- as.data.frame(metrics)
    colnames(metrics) <- c(
        "mean", "median", "G50 (assigned genes)",
        "G50 (all genes)", "O50 (assigned genes)",
        "O50 (all genes)", "All species OG", "sOGs"
    )
    GeneDiscoveRobject$runsData$overallMetrics <- metrics
    return(GeneDiscoveRobject)
}

#' Calculate statistics for selected files
#'
#' This function reads in a set of files that match the pattern "^N0" and calculates
#' various statistics based on the data in each file. The statistics include the median,
#' mean, variance, minimum, maximum, number of unique values, total number of rows, and
#' percentage of rows per unique value for a specific variable in each file.
#'
#' @param GeneDiscoveRobject An object of class GeneDiscoveR. Default is NULL.
#' @param cores The number of cores to use for parallel processing. Default is 1.
#' @return A GeneDiscoveR object with the calculated statistics added to the runsData$medians field.
#'
#' @import dplyr tidyr readr
#' @import doParallel
#' @import parallel
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
#' # Calculate median statistics
#' GeneDiscoveRobject <- calculate_median_statistics(GeneDiscoveRobject, cores = 8)
#' @export
calculate_median_statistics <- function(GeneDiscoveRobject = NULL, cores = 1) {
    if (is.null(GeneDiscoveRobject) || is.null(GeneDiscoveRobject$N0sDir)) {
        stop("Please provide the GeneDiscoveRobject parameter and make sure GeneDiscoveRobject$N0sDir is not null.")
    }
    medians <- NULL
    allMedians <- NULL

    # Parallelize the loop using 'foreach' and 'doParallel' packages
    if (cores > 1) {
        # Check if the number of threads is greater than the number of available cores
        if (cores > detectCores()) {
            stop("The number of cores cannot be greater than the number of available cores.")
        }

        # Create a cluster with the specified number of cores
        cl <- makeCluster(cores)
        registerDoParallel(cl)

        # Loop through the files in parallel
        medians <- foreach(file = list.files(path = GeneDiscoveRobject$N0sDir, pattern = "^N0", full.names = TRUE), .combine = rbind, .packages = c("dplyr", "readr")) %dopar% {
            table <- suppressMessages(read_tsv(file, progress = FALSE))
            N <- table %>%
                group_by(OG) %>%
                summarise(n = n())
            cbind(
                medianHOGs = median(N$n),
                meanHOGs = mean(N$n),
                varHOGs = var(N$n),
                minHOGs = min(N$n),
                maxHOGs = max(N$n),
                nOGstotal = length(unique(table$OG)),
                nHOGstotal = dim(table)[1],
                perPartitionByOG = ((dim(table)[1] * 100.0) / length(unique(table$OG))) - 100
            )
        }

        # Stop the cluster
        stopCluster(cl)
    } else {
        # Loop through the files sequentially
        for (file in list.files(path = GeneDiscoveRobject$N0sDir, pattern = "^N0", full.names = TRUE)) {
            table <- suppressMessages(read_tsv(file, progress = FALSE))
            N <- table %>%
                group_by(OG) %>%
                summarise(n = n())
            medians <- rbind(medians, cbind(
                medianHOGs = median(N$n),
                meanHOGs = mean(N$n),
                varHOGs = var(N$n),
                minHOGs = min(N$n),
                maxHOGs = max(N$n),
                nOGstotal = length(unique(table$OG)),
                nHOGstotal = dim(table)[1],
                perPartitionByOG = ((dim(table)[1] * 100.0) / length(unique(table$OG))) - 100
            ))
        }
    }

    medians <- medians %>%
        as.data.frame(medians)
    GeneDiscoveRobject$runsData$medians <- medians

    return(GeneDiscoveRobject)
}

#' Plot all species OGs and sOGs per inflation
#'
#' This function plots the all species OGs and sOGs per inflation.
#' If a run is active, a vertical dashed line is added to the plot on the active run's inflation value.
#' The function is only for executions with multiple inflation values.
#'
#' @param GeneDiscoveRobject An object of class "GeneDiscoveR" containing the necessary data for plotting.
#' @return An object of class "ggplot" representing the plot.
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
#' # Plot all species OGs and sOGs per inflation
#' plot_allSpeciesOGs_sOGs_per_inflation(GeneDiscoveRobject)
#'
#' # Plot all species OGs and sOGs per inflation with active run
#'
#' GeneDiscoveRobject <- set_run_active(GeneDiscoveRobject, InflationValue = 1.8)
#' plot_allSpeciesOGs_sOGs_per_inflation(GeneDiscoveRobject)
#'
#' @import ggplot2 ggsci
#' @export
plot_allSpeciesOGs_sOGs_per_inflation <- function(GeneDiscoveRobject = NULL) {
    if (is.null(GeneDiscoveRobject) || is.null(GeneDiscoveRobject$runsData$overallMetrics)) {
        stop("Please provide the GeneDiscoveRobject parameter and make sure GeneDiscoveRobject$runsData$overallMetrics are not null.")
    }
    if (length(GeneDiscoveRobject$Inflation) <= 1) {
        stop("This function is only for executions with multiple inflation values.")
    }
    .set_ggplot2_theme()
    plot <- ggplot(GeneDiscoveRobject$runsData$overallMetrics) +
        geom_line(aes(x = seq(GeneDiscoveRobject$Inflation[1], GeneDiscoveRobject$Inflation[2], GeneDiscoveRobject$Inflation[3]), y = as.numeric(`All species OG`), color = "All species OG"), alpha = 0.5) +
        geom_point(aes(x = seq(GeneDiscoveRobject$Inflation[1], GeneDiscoveRobject$Inflation[2], GeneDiscoveRobject$Inflation[3]), y = as.numeric(`All species OG`), color = "All species OG")) +
        geom_line(aes(x = seq(GeneDiscoveRobject$Inflation[1], GeneDiscoveRobject$Inflation[2], GeneDiscoveRobject$Inflation[3]), y = as.numeric(`sOGs`), color = "sOGs"), alpha = 0.5) +
        geom_point(aes(x = seq(GeneDiscoveRobject$Inflation[1], GeneDiscoveRobject$Inflation[2], GeneDiscoveRobject$Inflation[3]), y = as.numeric(`sOGs`), color = "sOGs")) +
        labs(title = "Metrics of Orthogroups", x = "Inflation value (I)", y = "Number of orthogroups (OG)") +
        scale_x_continuous(breaks = seq(GeneDiscoveRobject$Inflation[1], GeneDiscoveRobject$Inflation[2], GeneDiscoveRobject$Inflation[3]), labels = as.character(seq(GeneDiscoveRobject$Inflation[1], GeneDiscoveRobject$Inflation[2], GeneDiscoveRobject$Inflation[3]))) +
        scale_color_jama() +
        guides(color = guide_legend(title = "Metrics")) +
        scale_y_continuous(breaks = seq(0, max(GeneDiscoveRobject$runsData$overallMetrics$`All species OG`, GeneDiscoveRobject$runsData$overallMetrics$`sOGs`) + 10, 10), labels = as.character(seq(0, max(GeneDiscoveRobject$runsData$overallMetrics$`All species OG`, GeneDiscoveRobject$runsData$overallMetrics$`sOGs`) + 10, 10)))

    if (!is.null(GeneDiscoveRobject$RunActive)) {
        plot <- plot + geom_vline(xintercept = GeneDiscoveRobject$RunActive$InflationActive, linetype = "dashed", color = "coral")
    }

    return(plot)
}

#' Plot OGs and HOGs per inflation
#'
#' This function plots the number of orthogroups (OGs) and hierarchical orthogroups (HOGs) per inflation.
#' If a run is active, a vertical dashed line is added to the plot on the active run's inflation value.
#' The function is only for executions with multiple inflation values.
#'
#' @param GeneDiscoveRobject An object of class GeneDiscoveR. Default is NULL.
#' @return An object of class "ggplot" representing the plot.
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
#'
#' GeneDiscoveRobject <- calculate_median_statistics(GeneDiscoveRobject, cores = 8)
#'
#' # Plot number of orthogroups (OGs) and hierarchical orthogroups (HOGs) per inflation
#' plot_OGs_HOGs_per_inflation(GeneDiscoveRobject)
#'
#' # Plot number of orthogroups (OGs) and hierarchical orthogroups (HOGs) per inflation with active run
#'
#' GeneDiscoveRobject <- set_run_active(GeneDiscoveRobject, InflationValue = 1.8)
#' plot_OGs_HOGs_per_inflation(GeneDiscoveRobject)
#'
#' @import ggplot2 ggsci
#' @export
plot_OGs_HOGs_per_inflation <- function(GeneDiscoveRobject = NULL) {
    if (is.null(GeneDiscoveRobject) || is.null(GeneDiscoveRobject$runsData$medians)) {
        stop("Please provide the GeneDiscoveRobject parameter and make sure GeneDiscoveRobject$runsData$medians are not null.")
    }
    if (length(GeneDiscoveRobject$Inflation) <= 1) {
        stop("This function is only for executions with multiple inflation values.")
    }
    .set_ggplot2_theme()
    plot <- ggplot(GeneDiscoveRobject$runsData$medians) +
        geom_line(aes(x = seq(GeneDiscoveRobject$Inflation[1], GeneDiscoveRobject$Inflation[2], GeneDiscoveRobject$Inflation[3]), y = as.numeric(nOGstotal), color = "Number of OGs")) +
        geom_point(aes(x = seq(GeneDiscoveRobject$Inflation[1], GeneDiscoveRobject$Inflation[2], GeneDiscoveRobject$Inflation[3]), y = as.numeric(nOGstotal), color = "Number of OGs")) +
        geom_line(aes(x = seq(GeneDiscoveRobject$Inflation[1], GeneDiscoveRobject$Inflation[2], GeneDiscoveRobject$Inflation[3]), y = as.numeric(nHOGstotal), color = "Number of HOGs")) +
        geom_point(aes(x = seq(GeneDiscoveRobject$Inflation[1], GeneDiscoveRobject$Inflation[2], GeneDiscoveRobject$Inflation[3]), y = as.numeric(nHOGstotal), color = "Number of HOGs")) +
        labs(title = "Metrics of Orthogroups and HOGs", x = "Inflation value (I)", y = "Count") +
        scale_x_continuous(breaks = seq(GeneDiscoveRobject$Inflation[1], GeneDiscoveRobject$Inflation[2], GeneDiscoveRobject$Inflation[3]), labels = as.character(seq(GeneDiscoveRobject$Inflation[1], GeneDiscoveRobject$Inflation[2], GeneDiscoveRobject$Inflation[3]))) +
        scale_color_jama() +
        guides(color = guide_legend(title = "Metrics")) +
        scale_y_continuous(
            limits = c(as.numeric(min(GeneDiscoveRobject$runsData$medians$nOGstotal)), as.numeric(max(GeneDiscoveRobject$runsData$medians$nHOGstotal)) + 2000),
            breaks = seq(as.numeric(min(GeneDiscoveRobject$runsData$medians$nOGstotal)), as.numeric(max(GeneDiscoveRobject$runsData$medians$nHOGstotal)) + 2000, 2000),
            labels = as.character(seq(as.numeric(min(GeneDiscoveRobject$runsData$medians$nOGstotal)), as.numeric(max(GeneDiscoveRobject$runsData$medians$nHOGstotal)) + 2000, 2000))
        )

    if (!is.null(GeneDiscoveRobject$RunActive)) {
        plot <- plot + geom_vline(xintercept = GeneDiscoveRobject$RunActive$InflationActive, linetype = "dashed", color = "coral")
    }

    return(plot)
}
