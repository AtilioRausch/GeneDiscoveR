# Input files ---------------------------------------------------------------

# Run multiple findind of orthologs with different Inflation values
# Execute the following command in the terminal to run:
# snakemake --use-conda --cores 32 --snakefile run_orthofinder.py --configfile config.yaml

# Reorganize the data of the orthofinder results with the following bash scripts:
# bash/collect_rename_N0s_files.sh
# bash/collect_rename_overalls_files.sh

# Load the data into R -----------------------------------------------------

# Install and import PhenoR package
devtools::install_github("AtilioRausch/PhenoR")
library(PhenoR)
library(tidyverse)
library(ggsci)

# Directory where the data is located
overallsDir <- "/home/atilio/Escritorio/LiverwortGitHub/PhenoR/inst/extdata/Comparatives-1dot3-6/"
N0sDir <- "/home/atilio/Escritorio/LiverwortGitHub/PhenoR/inst/extdata/N0-1dot3-6/"
genomesTSV <- "inst/extdata/genomes.tsv"

PhenoRobject <- PhenoR(
    overallsDir = overallsDir,
    N0sDir = N0sDir,
    dataFile = genomesTSV,
    minInflation = 1.3,
    maxInflation = 6,
    stepInflation = 0.1
)

# system.file("extdata", "Comparatives-1.3-6", package = "PhenoR")
# system.file("extdata", "N0-1.3-6", package = "PhenoR")

# Calculate statistics for each execution of Orthofinder for selecting the best inflation value
PhenoRobject <- calculate_overall_statistics(PhenoRobject)
PhenoRobject <- calculate_median_statistics(PhenoRobject)

# Plot the data -----------------------------------------------------------
set_ggplot2_theme()
plotAllSpeciesOGs <- plot_allSpeciesOGs_sOGs_per_inflation(PhenoRobject)
plotOGsAndHOGs <- plot_OGs_HOGs_per_inflation(PhenoRobject)

# Read N0 file (result of Orthofinder)-----------------------------------------------------------
PhenoRobject <- set_run_active(PhenoRobject, InflationValue = 1.8)

PhenoRobject <- select_species_by_phenotype(
    PhenoRobject = PhenoRobject,
    columnPhenotype = `Oil-body-type`,
    columnOrthofinderID = `OrthofinderID`,
    type = "one in specialized cell"
)
PhenoRobject <- select_species_by_phenotype(
    PhenoRobject = PhenoRobject,
    columnPhenotype = `Oil-body-type`,
    columnOrthofinderID = `OrthofinderID`,
    type = "many in all cells"
)
PhenoRobject <- select_species_by_phenotype(
    PhenoRobject = PhenoRobject,
    columnPhenotype = `Oil-body-type`,
    columnOrthofinderID = `OrthofinderID`,
    type = "none"
)

# PhenoRobject <- clean_phenotypes(PhenoRobject)

PhenoRobject <- gene_identification_by_phenotype(formula = as.formula("one in specialized cell" ~ "many in all cells"), PhenoRobject = PhenoRobject, statistic = "Fisher", name = "PerType")
PhenoRobject <- select_genes_by_phenotype(PhenoRobject, thresholdPValue = 0.05, thresholdOddsRatio = 1)

filter0.05 <- countHOG %>% filter(pvalueFisher <= 0.05 && oddRatioFisher > 1)
F1anno <- map_annotation(tableOG_HOG = filter0.05, tableAnnotation = funcAnnotation)
F1anno <- split_annotation_per_gene(tableOG_HOG_Anno = F1anno)
F1anno <- filter_oneGene_per_OG_HOG(tableOG_HOG = F1anno)
filter0.05 <- merge(filter0.05, F1anno, by = c("OG", "HOG"), all = T)


g1 <- ggplot() +
    geom_point(aes(x = -log(countHOG$pvalueFisherPerType), y = countHOG$logOddRatioPerType), color = "#696D7D") +
    geom_point(aes(
        x = -log(countHOG$pvalueFisherPerType[which(str_detect(countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`, paste(markersClusterOBfiltavg_log2FC$gene, collapse = "|")))]),
        y = countHOG$logOddRatioPerType[which(str_detect(countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`, paste(markersClusterOBfiltavg_log2FC$gene, collapse = "|")))], color = "#7DCD85"
    )) +
    geom_point(aes(
        x = -log(countHOG$pvalueFisherPerType[which(str_detect(pattern = "Mp3g07510", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))]),
        y = countHOG$logOddRatioPerType[which(str_detect(pattern = "Mp3g07510", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))], color = "#E6AF2E", size = 1
    ), shape = 18) +
    geom_point(aes(
        x = -log(countHOG$pvalueFisherPerType[which(str_detect(pattern = "Mp6g08690", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))]),
        y = countHOG$logOddRatioPerType[which(str_detect(pattern = "Mp6g08690", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))], color = "#FF1B1C", size = 1
    ), shape = 18) +
    geom_point(aes(
        x = -log(countHOG$pvalueFisherPerType[which(str_detect(pattern = "Mp4g20670", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))]),
        y = countHOG$logOddRatioPerType[which(str_detect(pattern = "Mp4g20670", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))], color = "#BF4E30", size = 1
    ), shape = 18) +
    geom_point(aes(
        x = -log(countHOG$pvalueFisherPerType[which(str_detect(pattern = "Mp5g01530", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))]),
        y = countHOG$logOddRatioPerType[which(str_detect(pattern = "Mp5g01530", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))], size = 1, color = "#DCD6F7"
    ), shape = 18) +
    geom_point(aes(
        x = -log(countHOG$pvalueFisherPerType[which(str_detect(pattern = "Mp3g02320", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))]),
        y = countHOG$logOddRatioPerType[which(str_detect(pattern = "Mp3g02320", string = countHOG$`MpTAKv6-Marchantia_polymorpha_rudelaris`))], size = 1, color = "#F67E7D"
    ), shape = 18) +
    scale_color_manual("HOG of gene...",
        breaks = c("#7DCD85", "#E6AF2E", "#FF1B1C", "#BF4E30", "#DCD6F7", "#F67E7D"),
        values = c("#7DCD85", "#E6AF2E", "#FF1B1C", "#BF4E30", "#DCD6F7", "#F67E7D"),
        labels = c("In single cell", "MpMYB2", "MpERF13", "MpSYP12B", "MpNAC6", "MpC1HDZ")
    ) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dotted", col = "red") +
    geom_vline(xintercept = 2.9957, linetype = "dotted", col = "red") +
    annotate("text", x = 2.9957, y = 0, label = "p-value <= 0.05", vjust = 1.5, color = "black") +
    annotate("text", x = 8, y = 1, label = "OddRatio >= 1", color = "black") +
    ylab("") +
    xlab("-log(p-value)") +
    ggtitle("Fisher's Exact Test per Type of Oil-body")
