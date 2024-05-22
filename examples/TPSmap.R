# Input files ---------------------------------------------------------------

# Run multiple findind of orthologs with different Inflation values
# Execute the following command in the terminal to run:
# snakemake --use-conda --cores 32 --snakefile run_orthofinder.py --configfile config.yaml

# Reorganize the data of the orthofinder results with the following bash scripts:
# bash/collect_rename_N0s_files.sh
# bash/collect_rename_overalls_files.sh

# Load the data into R -----------------------------------------------------

# Install and import GeneDiscoveR package
devtools::install_github("AtilioRausch/GeneDiscoveR")
library(GeneDiscoveR)
library(tidyverse)
library(ggsci)
library(doParallel)

# Directory where the data is located
overallsDir <- "/home/atilio/Escritorio/LiverwortGitHub/GeneDiscoveR/inst/extdata/Comparatives-1dot3-6/"
N0sDir <- "/home/atilio/Escritorio/LiverwortGitHub/GeneDiscoveR/inst/extdata/N0-1dot3-6/"
genomesTSV <- "inst/extdata/genomes.tsv"

GeneDiscoveRobject <- GeneDiscoveR(
  overallsDir = overallsDir,
  N0sDir = N0sDir,
  dataFile = genomesTSV,
  minInflation = 1.3,
  maxInflation = 6,
  stepInflation = 0.1
)

# system.file("extdata", "Comparatives-1.3-6", package = "GeneDiscoveR")
# system.file("extdata", "N0-1.3-6", package = "GeneDiscoveR")

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
set_ggplot2_theme()
plotAllSpeciesOGs <- plot_allSpeciesOGs_sOGs_per_inflation(GeneDiscoveRobject)
plotOGsAndHOGs <- plot_OGs_HOGs_per_inflation(GeneDiscoveRobject)

# Read N0 file (result of Orthofinder)-----------------------------------------------------------
GeneDiscoveRobject <- set_run_active(GeneDiscoveRobject, InflationValue = 1.8)

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

# GeneDiscoveRobject <- select_species_by_phenotype(
#     GeneDiscoveRobject = GeneDiscoveRobject,
#     columnPhenotype = `Oil-body-type`,
#     columnOrthofinderID = `OrthofinderID`,
#     type = "noneOB"
# )

# GeneDiscoveRobject <- clean_phenotypes(GeneDiscoveRobject)library(shiny)
# To do
# GeneDiscoveRobject <- clean_identification(GeneDiscoveRobject)



overallsDir <- "/home/atilio/Escritorio/LiverwortGitHub/GeneDiscoveR/inst/extdata/Comparatives-1dot3-6/"
N0sDir <- "/home/atilio/Escritorio/LiverwortGitHub/GeneDiscoveR/inst/extdata/N0-1dot3-6/"
proteomesTSV <- "inst/extdata/genomes.tsv"
GeneDiscoveRobject <- GeneDiscoveR(
  overallsDir = overallsDir,
  N0sDir = N0sDir,
  dataFile = genomesTSV,
  minInflation = 1.3,
  maxInflation = 6,
  stepInflation = 0.1
)
GeneDiscoveRobject <- set_run_active(GeneDiscoveRobject, InflationValue = 1.8)
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
GeneDiscoveRobject <- select_species_by_phenotype(
  GeneDiscoveRobject = GeneDiscoveRobject,
  columnPhenotype = "Oil-body-type",
  columnID = "OrthofinderID",
  type = "noneOB"
)

t <- proc.time() # Inicia el cronómetro
GeneDiscoveRobject <- gene_identification_by_phenotype(
  GeneDiscoveRobject = GeneDiscoveRobject,
  formula = as.formula("noneOB ~ one_in_specialized_cell + many_in_all_cells"),
  statistic = "Fisher",
  name = "OBpresence",
  cores = 8
)
proc.time() - t
#    user  system elapsed                        1 core 4Ghz
#  61.421   0.003  61.432
#    user  system elapsed                        8 cores 4Ghz
# 107.199   4.617  16.453

# GeneDiscoveRobject <- gene_identification_by_phenotype(
#     GeneDiscoveRobject = GeneDiscoveRobject,
#     formula = as.formula("many_in_all_cells ~ one_in_specialized_cell"),
#     statistic = "Fisher",
#     name = "PerType"
# )

GeneDiscoveRobject <- select_genes_by_phenotype(GeneDiscoveRobject,
  pvalue = 0.05,
  oddsRatio = 1,
  sign = ">=",
  name = "OBpresence"
)

dim(GeneDiscoveRobject$FilteredGenes[[1]]$table)

View(GeneDiscoveRobject$RunActive$N0Active)
# GeneDiscoveRobject <- select_genes_by_phenotype(GeneDiscoveRobject,
#     pvalue = 0.05,
#     oddsRatio = 1,
#     sign = ">=",
#     name = "PerType"
# )

View(GeneDiscoveRobject$Identification)
View(GeneDiscoveRobject$FilteredGenes[[1]]$table)

GeneDiscoveRobject <- set_annotation_file(GeneDiscoveRobject, annotationFile = "inst/extdata/MpTak_v6.1_func_annotation_1line.tsv")

indexFilteredGenes <- select_filtered_gene_index(GeneDiscoveRobject, name = "OBpresence", pvalue = 0.05, oddsRatio = 1, sign = ">=")
GeneDiscoveRobject <- map_annotation(
  GeneDiscoveRobject = GeneDiscoveRobject,
  indexFilteredGenes = indexFilteredGenes,
  specieWithAnnotation = "MpTAKv6-Marchantia_polymorpha_rudelaris"
)

identification <- get_identification(GeneDiscoveRobject, name = "OBpresence")

filteredTable <- get_filtered_genes_table(GeneDiscoveRobject, name = "OBpresence", pvalue = 0.05, oddsRatio = 1, sign = ">=")
plot_GeneDiscoveR_volcano(GeneDiscoveRobject, name = "OBpresence")
run_web_app()



log(countHOG$oddsRatioFisherOBpresence)
View(GeneDiscoveRobject$FilteredGenes[[indexFilteredGenes]]$table)
countHOG <- GeneDiscoveRobject$RunActive$N0Active


library(plotly)
g1_plotly <- ggplotly(g1)
