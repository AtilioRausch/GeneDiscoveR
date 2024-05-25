# This script shows an example of how to use the GeneDiscoveR package with InParanoiDB-Phytozome data.
# For this analisys, you can download the data from the following link:
# https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=Phytozome
# From Phytozome, you can download the InParanoiDB data for the species you want to analyze.
# From Phytozome12/global_analysis/orthology/inParanoid/Phytozome: inParanoid_current.tar.gz
# Then, you can extract all files with you principal species vs the other species.


# Intall and load the GeneDiscoveR package
invisible(lapply(c("usethis", "devtools"), library, character.only = TRUE))
devtools::install_github("AtilioRausch/GeneDiscoveR")
library(GeneDiscoveR)

pairSpeciesDir <- system.file("extdata", "inparanoidbFiles", package = "GeneDiscoveR")
dataFile <- system.file("extdata", "inparanoidbFiles", "datafile.tsv", package = "GeneDiscoveR")
principalSpecie <- "Athalianacolumbia"
principalSpeciePrefix <- "AT"

# Create a GeneDiscoveR object with InParanoiDB-Phytozome data
GeneDiscoveRobject <- GeneDiscoveR(
    pairSpeciesDir = pairSpeciesDir,
    dataFile = dataFile,
    principalSpecie = principalSpecie,
    principalSpeciePrefix = "AT",
    orthologsTool = "InParanoiDB"
)

# Set active run with 8 cores. Here, we are using the InParanoiDB-Phytozome data.
# This function processes the data to form the necessary structures for the identification of genes.
GeneDiscoveRobject <- set_run_active(GeneDiscoveRobject, cores = 8)

# Select species by phenotype
GeneDiscoveRobject <- select_species_by_phenotype(
    GeneDiscoveRobject = GeneDiscoveRobject,
    columnPhenotype = "Example",
    columnID = "InparanoiDBID",
    type = "1"
)
GeneDiscoveRobject <- select_species_by_phenotype(
    GeneDiscoveRobject = GeneDiscoveRobject,
    columnPhenotype = "Example",
    columnID = "InparanoiDBID",
    type = "0"
)

# Identify genes by phenotype
GeneDiscoveRobject <- gene_identification_by_phenotype(
    GeneDiscoveRobject = GeneDiscoveRobject,
    formula = as.formula("0 ~ 1"),
    statistic = "Fisher",
    name = "Example",
    cores = 8
)

# Select genes by phenotype with p-value <= 0.1 and odds ratio > 1
GeneDiscoveRobject <- select_genes_by_phenotype(GeneDiscoveRobject,
    pvalue = 0.1,
    oddsRatio = 1,
    sign = ">",
    name = "Example"
)

# Set annotation file
annotationFile <- system.file("extdata", "Brassicaceae", "TAIR10_functional_descriptions", package = "GeneDiscoveR")
GeneDiscoveRobject <- set_annotation_file(GeneDiscoveRobject, annotationFile = annotationFile)
# Select filtered gene index
indexFilteredGenes <- select_filtered_gene_index(GeneDiscoveRobject, name = "Example", pvalue = 0.1, oddsRatio = 1, sign = ">")

# Map annotation
GeneDiscoveRobject <- map_annotation(
    GeneDiscoveRobject = GeneDiscoveRobject,
    indexFilteredGenes = indexFilteredGenes,
    specieWithAnnotation = "Athalianacolumbia",
    oneColumn = FALSE,
    sep = "\t"
)
# Visualize the filtered genes table with annotation
View(get_filtered_genes_table(GeneDiscoveRobject, name = "Example", pvalue = 0.1, oddsRatio = 1, sign = ">"))

# For exporting the filtered genes table with annotation
write_csv(get_filtered_genes_table(GeneDiscoveRobject, name = "Example", pvalue = 0.1, oddsRatio = 1, sign = ">"), "filtered_genes_table.csv")

# For exporting complete table with identification
write_csv(get_complete_table(GeneDiscoveRobject), "complete_table.csv")

# Plot the volcano plot
plot_genediscover_volcano(GeneDiscoveRobject, name = "Example")

# Run the GeneDiscoveR web app
run_genediscover_web_app()

# End of the script
