# Input files ---------------------------------------------------------------

# Run multiple findind of orthologs with different Inflation values
# Execute the following command in the terminal to run:
# snakemake --use-conda --cores 32 --snakefile run_orthofinder.py --configfile config.yaml

# Reorganize the data of the orthofinder results with the following bash scripts:
# bash/collect_rename_N0s_files.sh
# bash/collect_rename_overalls_files.sh

# Load the data into R -----------------------------------------------------

# Install and import GeneDiscoveR package
invisible(lapply(c("usethis", "devtools"), library, character.only = TRUE))
devtools::install_github("AtilioRausch/GeneDiscoveR", force = T)
library(GeneDiscoveR)

# Directory where the data is located
overallsDir <- system.file("extdata", "Comparatives-1dot3-6", package = "GeneDiscoveR")
N0sDir <- system.file("extdata", "N0-1dot3-6", package = "GeneDiscoveR")
dataFile <- system.file("extdata", "annotatedCDSs.tsv", package = "GeneDiscoveR")

GeneDiscoveRobject <- GeneDiscoveR(
  overallsDir = overallsDir,
  N0sDir = N0sDir,
  dataFile = dataFile,
  minInflation = 1.3,
  maxInflation = 6,
  stepInflation = 0.1,
  orthologsTool = "OrthoFinder"
)

# Calculate statistics for each execution of Orthofinder for selecting the best inflation value

GeneDiscoveRobject <- calculate_overall_statistics(GeneDiscoveRobject, cores = 8)
#    user  system elapsed                        1 core 4Ghz
#   1.409   0.028   1.423
#    user  system elapsed                        8 cores 4Ghz
#   0.089   0.018   1.579

GeneDiscoveRobject <- calculate_median_statistics(GeneDiscoveRobject, cores = 8)
#    user  system elapsed                        1 core 4Ghz
#  45.191   0.104  41.457
#    user  system elapsed                        8 cores 4Ghz
#   0.144   0.008  12.612

# Plot the data -----------------------------------------------------------
plotAllSpeciesOGs <- plot_allSpeciesOGs_sOGs_per_inflation(GeneDiscoveRobject)
plotOGsAndHOGs <- plot_OGs_HOGs_per_inflation(GeneDiscoveRobject)

plotAllSpeciesOGs
plotOGsAndHOGs
# Read N0 file (result of Orthofinder)-----------------------------------------------------------
GeneDiscoveRobject <- set_run_active(GeneDiscoveRobject, InflationValue = 1.8)

# Plot the data with the selected inflation value -----------------------------------------------------------
plotAllSpeciesOGs <- plot_allSpeciesOGs_sOGs_per_inflation(GeneDiscoveRobject)
plotOGsAndHOGs <- plot_OGs_HOGs_per_inflation(GeneDiscoveRobject)

plotAllSpeciesOGs
plotOGsAndHOGs

# Select species by phenotype -----------------------------------------------------------
GeneDiscoveRobject <- select_species_by_phenotype(
  GeneDiscoveRobject = GeneDiscoveRobject,
  columnPhenotype = "Oil-body-type",
  columnID = "OrthofinderID",
  type = "one_in_specialized_cell"
)
GeneDiscoveRobject <- select_species_by_phenotype(
  GeneDiscoveRobject = GeneDiscoveRobject,
  columnPhenotype = "Oil-body-type",
  columnID = "OrthofinderID",
  type = "many_in_all_cells"
)

# For clean phenotypes
# GeneDiscoveRobject <- clean_phenotypes(GeneDiscoveRobject)

# Identify genes by phenotype -----------------------------------------------------------
GeneDiscoveRobject <- gene_identification_by_phenotype(
  GeneDiscoveRobject = GeneDiscoveRobject,
  formula = as.formula("many_in_all_cells ~ one_in_specialized_cell"),
  statistic = "Fisher",
  name = "PerType",
  cores = 8
)
#    user  system elapsed                        1 core 4Ghz
#  61.421   0.003  61.432
#    user  system elapsed                        8 cores 4Ghz
# 107.199   4.617  16.453


# Select genes by phenotype with p-value <= 0.05 and odds ratio >= 1 -----------------------------------------------------------
GeneDiscoveRobject <- select_genes_by_phenotype(GeneDiscoveRobject,
  pvalue = 0.05,
  oddsRatio = 1,
  sign = ">=",
  name = "PerType"
)

# Map annotation-----------------------------------------------------------
# Import Marchantia polymorpha annotation file from MarpolBase v6.1
annotationFile <- system.file("extdata", "MpTak_v6.1_func_annotation_1line.tsv", package = "GeneDiscoveR")
GeneDiscoveRobject <- set_annotation_file(GeneDiscoveRobject, annotationFile = annotationFile)

# Select filtered gene index
indexFilteredGenes <- select_filtered_gene_index(GeneDiscoveRobject, name = "PerType", pvalue = 0.05, oddsRatio = 1, sign = ">=")

# Map annotation to the filtered genes. indexFilteredGenes is the index of the filtered genes, if NULL, the annotation is mapped to the complete table
GeneDiscoveRobject <- map_annotation(
  GeneDiscoveRobject = GeneDiscoveRobject,
  indexFilteredGenes = indexFilteredGenes,
  specieWithAnnotation = "MpTAKv6-Marchantia_polymorpha_rudelaris",
  oneColumn = TRUE
)
# Obtain the filtered genes table with annotation
filteredTable <- get_filtered_genes_table(GeneDiscoveRobject, name = "PerType", pvalue = 0.05, oddsRatio = 1, sign = ">=")

# Map annotation to the complete table. indexFilteredGenes is NULL
GeneDiscoveRobject <- map_annotation(
  GeneDiscoveRobject = GeneDiscoveRobject,
  indexFilteredGenes = NULL,
  specieWithAnnotation = "MpTAKv6-Marchantia_polymorpha_rudelaris",
  oneColumn = TRUE
)
# Obtain the complete table with annotation
completeTable <- get_complete_table(GeneDiscoveRobject)

# Plot the volcano plot -----------------------------------------------------------
volcano <- plot_genediscover_volcano(GeneDiscoveRobject, name = "PerType")

# Load the TPS genes -----------------------------------------------------------
TPSgenes <- read_tsv(system.file("extdata", "TPSgenes.tsv"), col_names = T)

# Diterpene
GeneID <- TPSgenes$GeneID[TPSgenes$TPStype == "Diterpene"]
OrthoFinderID <- TPSgenes$OrthofinderID[TPSgenes$TPStype == "Diterpene"]
Diterpene <- obtain_OG_from_gene(GeneDiscoveRobject, GeneID, OrthoFinderID)
Diterpene <- Diterpene %>% mutate(TPStype = "Diterpene")

# Bacterial
GeneID <- TPSgenes$GeneID[TPSgenes$TPStype == "Bacterial"]
OrthoFinderID <- TPSgenes$OrthofinderID[TPSgenes$TPStype == "Bacterial"]
Bacterial <- obtain_OG_from_gene(GeneDiscoveRobject, GeneID, OrthoFinderID)
Bacterial <- Bacterial %>% mutate(TPStype = "Bacterial")

# Fungi
GeneID <- TPSgenes$GeneID[TPSgenes$TPStype == "Fungi"]
OrthoFinderID <- TPSgenes$OrthofinderID[TPSgenes$TPStype == "Fungi"]
Fungi <- obtain_OG_from_gene(GeneDiscoveRobject, GeneID, OrthoFinderID)
Fungi <- Fungi %>% mutate(TPStype = "Fungi")

# Plot the volcano plot with the TPS genes -----------------------------------------------------------
TPScategories <- c("Diterpene", "Bacterial", "Fungi")
title <- "Fisher's Exact Test per Type of Oil-body with TPS genes"
plot_genediscover_detector_volcano(GeneDiscoveRobject,
  annotationTable = TPSs, title = title,
  name = "PerType", type = "TPStype",
  categories = TPScategories
)

# Run the GeneDiscoveR web app -----------------------------------------------------------
run_genediscover_web_app()

# End of the script -----------------------------------------------------------
