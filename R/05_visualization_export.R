run_web_app <- function() {
    launch_web_app()
}

get_identification <- function(GeneDiscoveRobject = NULL, name = NULL) {
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

plot_GeneDiscoveR_volcano <- function(GeneDiscoveRobject = NULL, name = NULL) {
    GeneDiscoveRidentification <- get_identification(GeneDiscoveRobject = GeneDiscoveRobject, name = name)
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
    g1 <- ggplot() +
        geom_point(aes(x = -log(table[[GeneDiscoveRidentification$columns[4]]]), y = table$logoddRatioFisher), color = "#696D7D") +
        coord_flip() +
        geom_hline(yintercept = 0, linetype = "dotted", col = "red") +
        geom_vline(xintercept = 2.9957, linetype = "dotted", col = "red") +
        annotate("text", x = 2.9957, y = 0, label = "p-value <= 0.05", vjust = 1.5, color = "black") +
        annotate("text", x = 8, y = 1, label = "OddRatio >= 1", color = "black") +
        ylab("") +
        xlab("-log(p-value)")
    return(g1)
}
# geom_point(aes(x = -log(countHOG$pvalueFisherOBpresence[which(str_detect(pattern = "Mp3g07510", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))]), y = log(countHOG$oddsRatioFisherOBpresence[which(str_detect(pattern = "Mp3g07510", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))]), color = "#E6AF2E", size = 1), shape = 18) +
#   geom_point(aes(x = -log(countHOG$pvalueFisherOBpresence[which(str_detect(pattern = "Mp6g08690", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))]), y = log(countHOG$oddsRatioFisherOBpresence[which(str_detect(pattern = "Mp6g08690", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))]), color = "#FF1B1C", size = 1), shape = 18) +
#   geom_point(aes(x = -log(countHOG$pvalueFisherOBpresence[which(str_detect(pattern = "Mp4g20670", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))]), y = log(countHOG$oddsRatioFisherOBpresence[which(str_detect(pattern = "Mp4g20670", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))]), color = "#BF4E30", size = 1), shape = 18) +
#   geom_point(aes(x = -log(countHOG$pvalueFisherOBpresence[which(str_detect(pattern = "Mp5g01530", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))]), y = log(countHOG$oddsRatioFisherOBpresence[which(str_detect(pattern = "Mp5g01530", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))]), size = 1, color = "#DCD6F7"), shape = 18) +
#   geom_point(aes(x = -log(countHOG$pvalueFisherOBpresence[which(str_detect(pattern = "Mp3g02320", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))]), y = log(countHOG$oddsRatioFisherOBpresence[which(str_detect(pattern = "Mp3g02320", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))]), size = 1, color = "#F67E7D"), shape = 18) +
#   scale_color_manual("HOG of gene...",
#     breaks = c("#7DCD85", "#E6AF2E", "#FF1B1C", "#BF4E30", "#DCD6F7", "#F67E7D"),
#     values = c("#7DCD85", "#E6AF2E", "#FF1B1C", "#BF4E30", "#DCD6F7", "#F67E7D"),
#     labels = c("In single cell", "MpMYB2", "MpERF13", "MpSYP12B", "MpNAC6", "MpC1HDZ")
#   )
get_overall_statistics <- function(GeneDiscoveRobject = NULL) {
    if (is.null(GeneDiscoveRobject)) {
        stop("Error: 'GeneDiscoveRobject' cannot be NULL.")
    }
    if (is.null(GeneDiscoveRobject$runsData$overallMetrics)) {
        stop("Error: 'GeneDiscoveRobject$runsData$overallMetrics' cannot be NULL.")
    }
    return(GeneDiscoveRobject$runsData$overallMetrics)
}

get_medians_statistics <- function(GeneDiscoveRobject = NULL) {
    if (is.null(GeneDiscoveRobject)) {
        stop("Error: 'GeneDiscoveRobject' cannot be NULL.")
    }
    if (is.null(GeneDiscoveRobject$runsData$overallMetrics)) {
        stop("Error: 'GeneDiscoveRobject$runsData$medians' cannot be NULL.")
    }
    return(GeneDiscoveRobject$runsData$medians)
}



