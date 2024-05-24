#' Calculate statistics and return metrics dataframe
#'
#' This function calculates statistics from multiple files and returns a dataframe.
#'
#' @param overallDirectory The directory where the statistics files are located.
#' @param Nrun The number of runs.
#' @return A dataframe containing the calculated statistics.
#' 
#' @examples
#' overallsDir <- system.file("extdata", "Comparatives-1.3-6"metrics")
#' metrics <- calculate_overall_statistics(overallDirectory = overallsDir, Nrun = 48)
#' 
#' @import dplyr
#' @import tidyr
#' @import readr
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
            cbind(stats[1, 2], stats[2, 2], stats[3, 2], stats[4, 2],
                  stats[5, 2], stats[6, 2], stats[7, 2], stats[8, 2])
        }
        
        # Stop the cluster
        stopCluster(cl)
    } else {
        # Loop through the files sequentially
        for (index in indexs) {
            stats <- paste0(GeneDiscoveRobject$overallsDir, "/Statistics_Overall_", index, ".tsv")
            stats <- suppressMessages(readr::read_tsv(stats, col_names = FALSE, skip = 10, n_max = 8))
            metrics <- rbind(metrics, cbind(stats[1, 2], stats[2, 2], stats[3, 2], stats[4, 2],
                                             stats[5, 2], stats[6, 2], stats[7, 2], stats[8, 2]))
        }
    }
    
    metrics <- as.data.frame(metrics)
    colnames(metrics) <- c("mean", "median", "G50 (assigned genes)",
                           "G50 (all genes)", "O50 (assigned genes)",
                           "O50 (all genes)", "All species OG", "sOGs")
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
#' @param nsDirectory The directory where the files are located. Default is NULL.
#' @return A data frame containing the calculated statistics for each file.
#'
#' @import dplyr
#' @import tidyr
#' @import readr
#'
#' @examples
#' N0sDir <- system.file("extdata", "N0-1.3-6"metrics")
#' medians <- calculate_median_statistics(nsDirectory = N0sDir)
#' 
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
            table <- suppressMessages(readr::read_tsv(file, progress = FALSE))
            N <- table %>% 
            dplyr::group_by(OG) %>% 
            dplyr::summarise(n = n())
            cbind(medianHOGs = median(N$n), 
                  meanHOGs = mean(N$n),
                  varHOGs = var(N$n),
                  minHOGs = min(N$n),
                  maxHOGs = max(N$n),
                  nOGstotal = length(unique(table$OG)), 
                  nHOGstotal = dim(table)[1],
                  perPartitionByOG = ((dim(table)[1]*100.0)/length(unique(table$OG)))-100)
        }
        
        # Stop the cluster
        stopCluster(cl)
    } else {
        # Loop through the files sequentially
        for (file in list.files(path = GeneDiscoveRobject$N0sDir, pattern = "^N0", full.names = TRUE)) {
            table <- suppressMessages(readr::read_tsv(file, progress = FALSE))
            N <- table %>% 
                dplyr::group_by(OG) %>% 
                dplyr::summarise(n = n())
            medians <- rbind(medians, cbind(medianHOGs = median(N$n), 
                                            meanHOGs = mean(N$n),
                                            varHOGs = var(N$n),
                                            minHOGs = min(N$n),
                                            maxHOGs = max(N$n),
                                            nOGstotal = length(unique(table$OG)), 
                                            nHOGstotal = dim(table)[1],
                                            perPartitionByOG = ((dim(table)[1]*100.0)/length(unique(table$OG)))-100))
        }
    }
    
    medians <- medians %>% 
        as.data.frame(medians)
    GeneDiscoveRobject$runsData$medians <- medians
    
    return(GeneDiscoveRobject)
}

#' Plot Metrics of Orthogroups
#'
#' This function plots the metrics of orthogroups using the provided dataframe.
#' @param metrics A dataframe containing the calculated statistics.
#' @param Nrun The number of runs. Default is 48.
#' @return An object of class "ggplot" representing the plot.
#' @examples
#' overallsDir <- system.file("extdata", "Comparatives-1.3-6"metrics")
#' metrics <- calculate_overall_statistics(overallDirectory = overallsDir, Nrun = 48)
#' plot <- plot_allSpeciesOGs_sOGs_per_inflation(metrics)
#' 
#' @import ggplot2
#' @export
plot_allSpeciesOGs_sOGs_per_inflation <- function(GeneDiscoveRobject=NULL) {
    if (is.null(GeneDiscoveRobject) || is.null(GeneDiscoveRobject$runsData$overallMetrics) || is.null(GeneDiscoveRobject$InflationLimits)) {
        stop("Please provide the GeneDiscoveRobject parameter and make sure GeneDiscoveRobject$runsData$overallMetrics and GeneDiscoveRobject$InflationLimits are not null.")
    }
    plot <- ggplot(GeneDiscoveRobject$runsData$overallMetrics) +
        geom_line(aes(x=seq(GeneDiscoveRobject$InflationLimits[1], GeneDiscoveRobject$InflationLimits[2], GeneDiscoveRobject$InflationLimits[3]), y=as.numeric(`All species OG`), color="All species OG"), alpha=0.5) +
        geom_point(aes(x=seq(GeneDiscoveRobject$InflationLimits[1], GeneDiscoveRobject$InflationLimits[2], GeneDiscoveRobject$InflationLimits[3]), y=as.numeric(`All species OG`), color="All species OG")) +
        geom_line(aes(x=seq(GeneDiscoveRobject$InflationLimits[1], GeneDiscoveRobject$InflationLimits[2], GeneDiscoveRobject$InflationLimits[3]), y=as.numeric(`sOGs`), color="sOGs"), alpha=0.5) +
        geom_point(aes(x=seq(GeneDiscoveRobject$InflationLimits[1], GeneDiscoveRobject$InflationLimits[2], GeneDiscoveRobject$InflationLimits[3]), y=as.numeric(`sOGs`), color="sOGs")) +
        labs(title="Metrics of Orthogroups", x="Inflation value (I)", y="Number of orthogroups (OG)") +
        scale_x_continuous(breaks=seq(GeneDiscoveRobject$InflationLimits[1], GeneDiscoveRobject$InflationLimits[2], GeneDiscoveRobject$InflationLimits[3]), labels=as.character(seq(GeneDiscoveRobject$InflationLimits[1], GeneDiscoveRobject$InflationLimits[2], GeneDiscoveRobject$InflationLimits[3]))) +
        scale_color_jama() +
        guides(color=guide_legend(title="Metrics")) +
        scale_y_continuous(breaks=seq(0, max(GeneDiscoveRobject$runsData$overallMetrics$`All species OG`, GeneDiscoveRobject$runsData$overallMetrics$`sOGs`) + 10, 10), labels=as.character(seq(0, max(GeneDiscoveRobject$runsData$overallMetrics$`All species OG`, GeneDiscoveRobject$runsData$overallMetrics$`sOGs`) + 10, 10)))
    
    if (!is.null(GeneDiscoveRobject$RunActive)) {
        plot <- plot + geom_vline(xintercept = GeneDiscoveRobject$RunActive$InflationActive, linetype = "dashed", color = "coral")
    }
    
    return(plot)
}

#' Plot Metrics of Orthogroups and HOGs
#'
#' This function plots the metrics of orthogroups and HOGs using the provided dataframe.
#'
#' @param medians A dataframe containing the calculated statistics.
#' @param Nrun The number of runs.
#' @return An object of class "ggplot" representing the plot.
#' @examples
#' N0sDir <- system.file("extdata", "N0-1.3-6"metrics")
#' medians <- calculate_median_statistics(nsDirectory = N0sDir)
#' plot <- plot_OGs_HOGs_per_inflation(medians, Nrun = 48)
#'
#' @import ggplot2
#' @import ggsci
#' @export
plot_OGs_HOGs_per_inflation <- function(GeneDiscoveRobject=NULL) {
    if (is.null(GeneDiscoveRobject) || is.null(GeneDiscoveRobject$runsData$medians) || is.null(GeneDiscoveRobject$InflationLimits)) {
        stop("Please provide the GeneDiscoveRobject parameter and make sure GeneDiscoveRobject$runsData$medians and GeneDiscoveRobject$InflationLimits are not null.")
    }
    plot <- ggplot(GeneDiscoveRobject$runsData$medians) +
        geom_line(aes(x = seq(GeneDiscoveRobject$InflationLimits[1], GeneDiscoveRobject$InflationLimits[2], GeneDiscoveRobject$InflationLimits[3]), y = as.numeric(nOGstotal), color = "Number of OGs")) +
        geom_point(aes(x = seq(GeneDiscoveRobject$InflationLimits[1], GeneDiscoveRobject$InflationLimits[2], GeneDiscoveRobject$InflationLimits[3]), y = as.numeric(nOGstotal), color = "Number of OGs")) +
        geom_line(aes(x = seq(GeneDiscoveRobject$InflationLimits[1], GeneDiscoveRobject$InflationLimits[2], GeneDiscoveRobject$InflationLimits[3]), y = as.numeric(nHOGstotal), color = "Number of HOGs")) +
        geom_point(aes(x = seq(GeneDiscoveRobject$InflationLimits[1], GeneDiscoveRobject$InflationLimits[2], GeneDiscoveRobject$InflationLimits[3]), y = as.numeric(nHOGstotal), color = "Number of HOGs")) +
        labs(title = "Metrics of Orthogroups and HOGs", x = "Inflation value (I)", y = "Count") +
        scale_x_continuous(breaks = seq(GeneDiscoveRobject$InflationLimits[1], GeneDiscoveRobject$InflationLimits[2], GeneDiscoveRobject$InflationLimits[3]), labels = as.character(seq(GeneDiscoveRobject$InflationLimits[1], GeneDiscoveRobject$InflationLimits[2], GeneDiscoveRobject$InflationLimits[3]))) +
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
#'
#' @examples
#' set_ggplot2_theme()
#'
#' @export
set_ggplot2_theme <- function() {
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
        axis.title.y = element_text(vjust = 1,size = 15, angle = 90),
        legend.key = element_blank(),
        panel.border = element_rect(
            colour = "#77767b",
            fill = NA,
            linewidth = 0.5
        )
    )
    theme_set(theme)
}















