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

# Directory where the data is located
overallsDir <- system.file("extdata", "Comparatives-1.3-6", package = "PhenoR")
N0sDir <- system.file("extdata", "N0-1.3-6")

# Calculate statistics for each execution of Orthofinder for selecting the best inflation value
metrics <- PhenoR::calculate_overall_statistics(overallDirectory = overallsDir, Nrun = 48)
medians <- PhenoR::calculate_median_statistics(nsDirectory = N0sDir)

# Plot the data -----------------------------------------------------------
PhenoR::set_ggplot2_theme()
allSpeciesOGs <- PhenoR::plot_allSpeciesOGs_sOGs_per_inflation(metrics)
OGsAndHOGs <- PhenoR::plot_OGs_HOGs_per_inflation(medians, Nrun = 48)

# Read N0 file (result of Orthofinder)-----------------------------------------------------------

N0 <- read_tsv(file = "N0.tsv")
genomesTSV <- refactor_genomes_table(file = "/home/atilio/Escritorio/Code/genomes.tsv")
oneInSpecializedCellG <- select_species_by_OBtype(genomes = genomesTSV, type = "one in specialized cell", origin = "G")
manyInAllCellsTG <- select_species_by_OBtype(genomes = genomesTSV, type = "many in all cells", origin = "G")
onePerCellG <- select_species_by_OBtype(genomes = genomesTSV, type = "one per cell", origin = "G")
noneG <- select_species_by_OBtype(genomes = genomesTSV, type = "none", origin = "G")

countHOG <- countHOG %>%
    rowwise() %>%
    mutate(
        NInSpecializedCell = sum(!is.na(c_across(all_of(oneInSpecializedCellG)))),
        NmanyInAllCells = sum(!is.na(c_across(all_of(manyInAllCellsTG))))
    )
countHOG <- countHOG %>% mutate(
    pvalueFisherPerType = fisher.test(matrix(c(NInSpecializedCell, NmanyInAllCells, length(oneInSpecializedCellG) - NInSpecializedCell, length(manyInAllCellsTG) - NmanyInAllCells), nrow = 2, byrow = T))$p.value,
    oddRatioFisherPerType = fisher.test(matrix(c(NInSpecializedCell, NmanyInAllCells, length(oneInSpecializedCellG) - NInSpecializedCell, length(manyInAllCellsTG) - NmanyInAllCells), nrow = 2, byrow = T))$estimate
)
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
