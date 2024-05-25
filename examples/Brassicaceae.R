# This script is an example of how to use the GeneDiscoveR package to identify genes associated with a phenotype in the Brassicaceae family.

# For generate the input data for this example, you can use the following commands:
# Install OrthoFinder
# conda install -c bioconda orthofinder
# One execution of OrthoFinder with default parameters
# orthofinder -f OrthoFinder/ExampleData
# In the ExampleData directory, you locate the annotated coding sequences of your species for analysis.

# Import and install the GeneDiscoveR package
invisible(lapply(c("usethis", "devtools"), library, character.only = TRUE))
devtools::install_github("AtilioRausch/GeneDiscoveR")
library(GeneDiscoveR)

N0Dir <- system.file("extdata", "Brassicaceae", package = "GeneDiscoveR")
dataTSV <- system.file("extdata", "Brassicaceae", "table_traits_selfcomp.tsv", package = "GeneDiscoveR")

# Create a GeneDiscoveR object with Brassicaceae data from one execution of OrthoFinder
# In this case, because we are unique execution, we use the uniqueInflation parameter and GeneDiscoveR set automatically the active run.
GeneDiscoveRobject <- GeneDiscoveR(
    N0sDir = N0Dir,
    dataFile = dataTSV,
    uniqueInflation = 1.5,
    orthologsTool = "OrthoFinder"
)

# Select species by phenotype. You can perform this step with different phenotypes.
GeneDiscoveRobject <- select_species_by_phenotype(
    GeneDiscoveRobject = GeneDiscoveRobject,
    columnPhenotype = "Self-compatible",
    columnID = "OrthofinderID",
    type = "0"
)
GeneDiscoveRobject <- select_species_by_phenotype(
    GeneDiscoveRobject = GeneDiscoveRobject,
    columnPhenotype = "Self-compatible",
    columnID = "OrthofinderID",
    type = "1"
)

# Identify genes by phenotype. You can perform this step with different phenotypes.
GeneDiscoveRobject <- gene_identification_by_phenotype(
    GeneDiscoveRobject = GeneDiscoveRobject,
    formula = as.formula("0 ~ 1"),
    statistic = "Fisher",
    name = "Self-compatible",
    cores = 8
)

# Select genes by phenotype with p-value <= 0.05 and odds ratio >= 1
GeneDiscoveRobject <- select_genes_by_phenotype(GeneDiscoveRobject,
    pvalue = 0.05,
    oddsRatio = 1,
    sign = ">=",
    name = "Self-compatible"
)

# Set annotation file
# Import Arabidopsis thaliana annotation file from TAIR10
annotationFile <- system.file("extdata", "Brassicaceae", "TAIR10_functional_descriptions", package = "GeneDiscoveR")
GeneDiscoveRobject <- set_annotation_file(GeneDiscoveRobject, annotationFile = annotationFile)
indexFilteredGenes <- select_filtered_gene_index(GeneDiscoveRobject, name = "Self-compatible", pvalue = 0.05, oddsRatio = 1, sign = ">=")

# Map annotation to the filtered genes. indexFilteredGenes is the index of the filtered genes, if NULL, the annotation is mapped to the complete table
GeneDiscoveRobject <- map_annotation(
    GeneDiscoveRobject = GeneDiscoveRobject,
    indexFilteredGenes = indexFilteredGenes,
    specieWithAnnotation = "Athaliana_447_Araport11.protein_primaryTranscriptOnly",
    oneColumn = FALSE
)
# Map annotation to the complete table. indexFilteredGenes is NULL
GeneDiscoveRobject <- map_annotation(
    GeneDiscoveRobject = GeneDiscoveRobject,
    indexFilteredGenes = NULL,
    specieWithAnnotation = "Athaliana_447_Araport11.protein_primaryTranscriptOnly",
    oneColumn = FALSE
)

# Show the filtered genes table
View(get_filtered_genes_table(GeneDiscoveRobject, name = "Self-compatible", pvalue = 0.05, oddsRatio = 1, sign = ">="))

# plot volcano
plot_genediscover_volcano(GeneDiscoveRobject, name = "Self-compatible")

# Export the filtered genes table with annotation
write_csv(get_filtered_genes_table(GeneDiscoveRobject, name = "Self-compatible", pvalue = 0.05, oddsRatio = 1, sign = ">="), "/home/atilio/Descargas/Results_Jun02/filtered_genes_table-2.csv")

# Export the complete table with identification
write_csv(get_complete_table(GeneDiscoveRobject), "/home/path/to/complete_table-2.csv")

# Run the GeneDiscoveR web app
run_genediscover_web_app()

# End of the script
