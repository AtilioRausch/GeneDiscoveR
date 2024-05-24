devtools::install_github("AtilioRausch/GeneDiscoveR")
library(GeneDiscoveR)
library(tidyverse)
library(ggsci)

# Directory where the data is located
N0Dir <- "/home/atilio/Escritorio/Results_Jun02/"
dataTSV <- "/home/atilio/Escritorio/Results_Jun02/table_traits_selfcomp.tsv"

GeneDiscoveRobject <- GeneDiscoveR(
    N0sDir = N0Dir,
    dataFile = dataTSV,
    uniqueInflation = 1.5
)

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

GeneDiscoveRobject <- gene_identification_by_phenotype(
    GeneDiscoveRobject = GeneDiscoveRobject,
    formula = as.formula("0 ~ 1"),
    statistic = "Fisher",
    name = "Self-compatible",
    cores = 8
)

GeneDiscoveRobject <- select_genes_by_phenotype(GeneDiscoveRobject,
    pvalue = 0.05,
    oddsRatio = 1,
    sign = ">=",
    name = "Self-compatible"
)

GeneDiscoveRobject <- set_annotation_file(GeneDiscoveRobject, annotationFile = "/home/atilio/Descargas/Results_Jun02/TAIR10_functional_descriptions")
indexFilteredGenes <- select_filtered_gene_index(GeneDiscoveRobject, name = "Self-compatible", pvalue = 0.05, oddsRatio = 1, sign = ">=")
GeneDiscoveRobject <- map_annotation(
    GeneDiscoveRobject = GeneDiscoveRobject,
    indexFilteredGenes = indexFilteredGenes,
    specieWithAnnotation = "Athaliana_447_Araport11.protein_primaryTranscriptOnly",
    oneColumn = FALSE
)
View(get_filtered_genes_table(GeneDiscoveRobject, name = "Self-compatible", pvalue = 0.05, oddsRatio = 1, sign = ">="))

filteredTable <- get_filtered_genes_table(GeneDiscoveRobject, name = "Self-compatible", pvalue = 0.05, oddsRatio = 1, sign = ">=")
write_csv(filteredTable, "/home/atilio/Descargas/Results_Jun02/filteredTable-2.csv")

run_genediscover_web_app()
