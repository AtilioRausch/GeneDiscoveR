devtools::install_github("AtilioRausch/PhenoR")
library(PhenoR)
library(tidyverse)
library(ggsci)

# Directory where the data is located
N0Dir <- "/home/atilio/Descargas/Results_Jun02/"
dataTSV <- "/home/atilio/Descargas/Results_Jun02/table_traits_selfcomp.tsv"

PhenoRobject <- PhenoR(
    N0sDir = N0Dir,
    dataFile = dataTSV,
    uniqueInflation = 1.5
)

PhenoRobject <- select_species_by_phenotype(
    PhenoRobject = PhenoRobject,
    columnPhenotype = "Self-compatible",
    columnID = "OrthofinderID",
    type = "0"
)
PhenoRobject <- select_species_by_phenotype(
    PhenoRobject = PhenoRobject,
    columnPhenotype = "Self-compatible",
    columnID = "OrthofinderID",
    type = "1"
)

PhenoRobject <- gene_identification_by_phenotype(
    PhenoRobject = PhenoRobject,
    formula = as.formula("0 ~ 1"),
    statistic = "Fisher",
    name = "Self-compatible",
    cores = 8
)

PhenoRobject <- select_genes_by_phenotype(PhenoRobject,
    pvalue = 0.05,
    oddsRatio = 1,
    sign = ">=",
    name = "Self-compatible"
)

PhenoRobject <- set_annotation_file(PhenoRobject, annotationFile = "/home/atilio/Descargas/Results_Jun02/TAIR10_functional_descriptions")
indexFilteredGenes <- select_filtered_gene_index(PhenoRobject, name = "Self-compatible", pvalue = 0.05, oddsRatio = 1, sign = ">=")
PhenoRobject <- map_annotation(
    PhenoRobject = PhenoRobject,
    indexFilteredGenes = indexFilteredGenes,
    specieWithAnnotation = "Athaliana_447_Araport11.protein_primaryTranscriptOnly",
    oneColumn = FALSE
)
View(get_filtered_genes_table(PhenoRobject, name = "Self-compatible", pvalue = 0.05, oddsRatio = 1, sign = ">="))

filteredTable <- get_filtered_genes_table(PhenoRobject, name = "Self-compatible", pvalue = 0.05, oddsRatio = 1, sign = ">=")
write_csv(filteredTable, "/home/atilio/Descargas/Results_Jun02/filteredTable.csv")


run_web_app()
￼
