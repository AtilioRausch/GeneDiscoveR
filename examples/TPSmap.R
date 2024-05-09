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
allSpeciesOGs <- PhenoR::plot_allSpeciesOGs_sOGs_per_inflation(metrics)
OGsAndHOGs <- PhenoR::plot_OGs_HOGs_per_inflation(medians, Nrun = 48)
