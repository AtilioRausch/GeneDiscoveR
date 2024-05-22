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
library(doParallel)

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

PhenoRobject <- calculate_overall_statistics(PhenoRobject, cores = 8)
#    user  system elapsed                        1 core 4Ghz
#   1.409   0.028   1.423
#    user  system elapsed                        8 cores 4Ghz
#   0.089   0.018   1.579

PhenoRobject <- calculate_median_statistics(PhenoRobject, cores = 8)
#    user  system elapsed                        1 core 4Ghz
#  45.191   0.104  41.457
#    user  system elapsed                        8 cores 4Ghz
#   0.144   0.008  12.612

# Plot the data -----------------------------------------------------------
set_ggplot2_theme()
plotAllSpeciesOGs <- plot_allSpeciesOGs_sOGs_per_inflation(PhenoRobject)
plotOGsAndHOGs <- plot_OGs_HOGs_per_inflation(PhenoRobject)

# Read N0 file (result of Orthofinder)-----------------------------------------------------------
PhenoRobject <- set_run_active(PhenoRobject, InflationValue = 1.8)

PhenoRobject <- select_species_by_phenotype(
  PhenoRobject = PhenoRobject,
  columnPhenotype = "Oil-body-type",
  columnID = "OrthofinderID",
  type = "one_in_specialized_cell"
)
PhenoRobject <- select_species_by_phenotype(
  PhenoRobject = PhenoRobject,
  columnPhenotype = "Oil-body-type",
  columnID = "OrthofinderID",
  type = "many_in_all_cells"
)

# PhenoRobject <- select_species_by_phenotype(
#     PhenoRobject = PhenoRobject,
#     columnPhenotype = `Oil-body-type`,
#     columnOrthofinderID = `OrthofinderID`,
#     type = "noneOB"
# )

# PhenoRobject <- clean_phenotypes(PhenoRobject)library(shiny)
# To do
# PhenoRobject <- clean_identification(PhenoRobject)



overallsDir <- "/home/atilio/Escritorio/LiverwortGitHub/PhenoR/inst/extdata/Comparatives-1dot3-6/"
N0sDir <- "/home/atilio/Escritorio/LiverwortGitHub/PhenoR/inst/extdata/N0-1dot3-6/"
proteomesTSV <- "inst/extdata/genomes.tsv"
PhenoRobject <- PhenoR(
  overallsDir = overallsDir,
  N0sDir = N0sDir,
  dataFile = genomesTSV,
  minInflation = 1.3,
  maxInflation = 6,
  stepInflation = 0.1
)
PhenoRobject <- set_run_active(PhenoRobject, InflationValue = 1.8)
PhenoRobject <- select_species_by_phenotype(
  PhenoRobject = PhenoRobject,
  columnPhenotype = "Oil-body-type",
  columnID = "OrthofinderID",
  type = "one_in_specialized_cell"
)
PhenoRobject <- select_species_by_phenotype(
  PhenoRobject = PhenoRobject,
  columnPhenotype = "Oil-body-type",
  columnID = "OrthofinderID",
  type = "many_in_all_cells"
)
PhenoRobject <- select_species_by_phenotype(
  PhenoRobject = PhenoRobject,
  columnPhenotype = "Oil-body-type",
  columnID = "OrthofinderID",
  type = "noneOB"
)

t <- proc.time() # Inicia el cronómetro
PhenoRobject <- gene_identification_by_phenotype(
  PhenoRobject = PhenoRobject,
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

# PhenoRobject <- gene_identification_by_phenotype(
#     PhenoRobject = PhenoRobject,
#     formula = as.formula("many_in_all_cells ~ one_in_specialized_cell"),
#     statistic = "Fisher",
#     name = "PerType"
# )

PhenoRobject <- select_genes_by_phenotype(PhenoRobject,
  pvalue = 0.05,
  oddsRatio = 1,
  sign = ">=",
  name = "OBpresence"
)

dim(PhenoRobject$FilteredGenes[[1]]$table)

View(PhenoRobject$RunActive$N0Active)
# PhenoRobject <- select_genes_by_phenotype(PhenoRobject,
#     pvalue = 0.05,
#     oddsRatio = 1,
#     sign = ">=",
#     name = "PerType"
# )

View(PhenoRobject$Identification)
View(PhenoRobject$FilteredGenes[[1]]$table)

PhenoRobject <- set_annotation_file(PhenoRobject, annotationFile = "inst/extdata/MpTak_v6.1_func_annotation_1line.tsv")

indexFilteredGenes <- select_filtered_gene_index(PhenoRobject, name = "OBpresence", pvalue = 0.05, oddsRatio = 1, sign = ">=")
PhenoRobject <- map_annotation(
  PhenoRobject = PhenoRobject,
  indexFilteredGenes = indexFilteredGenes,
  specieWithAnnotation = "MpTAKv6-Marchantia_polymorpha_rudelaris"
)

identification <- get_identification(PhenoRobject, name = "OBpresence")

filteredTable <- get_filtered_genes_table(PhenoRobject, name = "OBpresence", pvalue = 0.05, oddsRatio = 1, sign = ">=")
plot_phenor_volcano(PhenoRobject, name = "OBpresence")
run_web_app()



log(countHOG$oddsRatioFisherOBpresence)
View(PhenoRobject$FilteredGenes[[indexFilteredGenes]]$table)
countHOG <- PhenoRobject$RunActive$N0Active


library(plotly)
g1_plotly <- ggplotly(g1)
