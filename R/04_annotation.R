#' Set the annotation file for a GeneDiscoveR object
#'
#' This function sets the annotation file for a GeneDiscoveR object.
#' The annotation file contains information about the genes, such as their names, IDs, and functional annotations.
#' The annotation file should be in tab-separated values (TSV) format, with the first column containing the gene IDs and the subsequent columns containing the annotations in one column or multiple columns.
#'
#' @param GeneDiscoveRobject A GeneDiscoveR object to which the annotation file will be set.
#' @param annotationFile The path to the annotation file.
#'
#' @return The updated GeneDiscoveR object with the annotation file set.
#' @export
#'
#' @examples
#' # Create a GeneDiscoveR object
#' N0sDir <- system.file("extdata", "N0-1dot3-6", package = "GeneDiscoveR")
#' overallsDir <- system.file("extdata", "Comparatives-1dot3-6", package = "GeneDiscoveR")
#' dataFile <- system.file("extdata", "annotatedCDSs.tsv", package = "GeneDiscoveR")
#' minInflation <- 1.3
#' maxInflation <- 6
#' stepInflation <- 0.1
#'
#' GeneDiscoveRobject <- GeneDiscoveR(overallsDir = overallsDir, N0sDir = N0sDir, dataFile = dataFile, minInflation = minInflation, maxInflation = maxInflation, stepInflation = stepInflation)
#'
#' # Set the annotation file
#' annotationFile <- system.file("extdata", "MpTak_v6.1_func_annotation_1line.tsv", package = "GeneDiscoveR")
#' GeneDiscoveRobject <- set_annotation_file(GeneDiscoveRobject, annotationFile = annotationFile)
set_annotation_file <- function(GeneDiscoveRobject = NULL, annotationFile = NULL) {
    if (is.null(GeneDiscoveRobject) || is.null(annotationFile)) {
        stop("Error: 'GeneDiscoveRobject' and 'annotationFile' cannot be NULL.")
    }
    GeneDiscoveRobject$AnnotationFile <- annotationFile
    return(GeneDiscoveRobject)
}
#' Map gene annotation
#'
#' This function maps gene annotation to the filtered genes in a GeneDiscoveR object.
#' The annotation file should contain in the first column the gene IDs and in the subsequent columns the annotations in one column or multiple columns.
#' If the annotation is in one column, the function will split the annotation into multiple columns by \code{sep} separator.
#' If the annotation is in multiple columns, the function will assign the annotation to the filtered genes based on the gene IDs.
#'
#' @param GeneDiscoveRobject A GeneDiscoveR object.
#' @param indexFilteredGenes The index of the filtered genes in the GeneDiscoveR object. If NULL, the annotation is mapped to the complete table.
#' @param oneColumn Logical value indicating whether to split the annotation into multiple columns (default is TRUE).
#' @param sep The separator to use when splitting the annotation into multiple columns (default is ";").
#' @param specieWithAnnotation The species with the annotation to be mapped (default is "MpTAKv6-Marchantia_polymorpha_rudelaris").
#'
#' @return The updated GeneDiscoveR object with the mapped gene annotation.
#'
#' @examples
#' # Create a GeneDiscoveR object
#' N0sDir <- system.file("extdata", "N0-1dot3-6", package = "GeneDiscoveR")
#' overallsDir <- system.file("extdata", "Comparatives-1dot3-6", package = "GeneDiscoveR")
#' dataFile <- system.file("extdata", "annotatedCDSs.tsv", package = "GeneDiscoveR")
#' minInflation <- 1.3
#' maxInflation <- 6
#' stepInflation <- 0.1
#'
#' GeneDiscoveRobject <- GeneDiscoveR(overallsDir = overallsDir, N0sDir = N0sDir, dataFile = dataFile, minInflation = minInflation, maxInflation = maxInflation, stepInflation = stepInflation)
#'
#' # Set active run
#' GeneDiscoveRobject <- set_run_active(GeneDiscoveRobject, InflationValue = 1.8)
#'
#' # Select species by phenotype
#' GeneDiscoveRobject <- select_species_by_phenotype(GeneDiscoveRobject = GeneDiscoveRobject, columnPhenotype = "Oil-body-type", columnID = "OrthofinderID", type = "one_in_specialized_cell")
#' GeneDiscoveRobject <- select_species_by_phenotype(GeneDiscoveRobject = GeneDiscoveRobject, columnPhenotype = "Oil-body-type", columnID = "OrthofinderID", type = "many_in_all_cells")
#'
#' # Gene identification by phenotype
#' GeneDiscoveRobject <- gene_identification_by_phenotype(GeneDiscoveRobject = GeneDiscoveRobject, formula = as.formula("one_in_specialized_cell ~ many_in_all_cells"), statistic = "Fisher", name = "PerOBtype", cores = 8)
#'
#' # Select genes by phenotype
#' GeneDiscoveRobject <- select_genes_by_phenotype(GeneDiscoveRobject, pvalue = 0.05, oddsRatio = 1, sign = ">=", name = "PerOBtype")
#'
#' # Annotation file with annotations in a one-line format separated by ";"
#' annotationFile <- system.file("extdata", "MpTak_v6.1_func_annotation_1line.tsv", package = "GeneDiscoveR")
#' # Set the annotation file
#' GeneDiscoveRobject <- set_annotation_file(GeneDiscoveRobject, annotationFile = annotationFile)
#' # Obtain the index of the filtered genes
#' indexFilteredGenes <- select_filtered_gene_index(GeneDiscoveRobject, name = "PerOBtype", pvalue = 0.05, oddsRatio = 1, sign = ">=")
#'
#' # Map gene annotation to the filtered genes
#' GeneDiscoveRobject <- map_annotation(GeneDiscoveRobject = GeneDiscoveRobject, indexFilteredGenes = indexFilteredGenes, specieWithAnnotation = "MpTAKv6-Marchantia_polymorpha_rudelaris", oneColumn = TRUE, sep = ";")
#' # Output: A GeneDiscoveR object with the mapped gene annotation
#'
#' # Map gene annotation to the complete table
#' GeneDiscoveRobject <- map_annotation(GeneDiscoveRobject = GeneDiscoveRobject, specieWithAnnotation = "MpTAKv6-Marchantia_polymorpha_rudelaris", oneColumn = TRUE, sep = ";")
#' # Output: A GeneDiscoveR object with the mapped gene annotation
#'
#' # Obtain the filtered genes table with annotation
#' filteredTable <- get_filtered_genes_table(GeneDiscoveRobject, name = "PerOBtype", pvalue = 0.05, oddsRatio = 1, sign = ">=")
#' # Obtain the complete table with annotation
#' completeTable <- get_complete_table(GeneDiscoveRobject)
#' # View the filtered genes table with the mapped gene annotation
#' filteredGenes <- get_filtered_genes_table(GeneDiscoveRobject, name = "PerOBtype", pvalue = 0.05, oddsRatio = 1, sign = ">=")
#' # Output: A table with the filtered genes and their annotations
#' @export
map_annotation <- function(GeneDiscoveRobject = NULL, indexFilteredGenes = NULL, oneColumn = TRUE, sep = ";", specieWithAnnotation = "MpTAKv6-Marchantia_polymorpha_rudelaris") {
    annotation <- .map_gene_annotation(GeneDiscoveRobject = GeneDiscoveRobject, indexFilteredGenes = indexFilteredGenes, specieWithAnnotation = specieWithAnnotation)
    if (oneColumn) {
        annotation <- .split_annotation_per_gene(tableOG_HOG_Anno = annotation, sep = sep)
    }
    if (GeneDiscoveRobject$OrthologsTool == "OrthoFinder") {
        annotation <- .filter_oneGene_per_OG_HOG(tableOG_HOG = annotation)
    }
    if (GeneDiscoveRobject$OrthologsTool == "InParanoiDB") {
        if (is.null(indexFilteredGenes)) {
            GeneDiscoveRobject$RunActive$N0Active <- merge(GeneDiscoveRobject$RunActive$N0Active, annotation, by = "OG", all = T)
        } else {
            GeneDiscoveRobject$FilteredGenes[[indexFilteredGenes]]$table <- merge(GeneDiscoveRobject$FilteredGenes[[indexFilteredGenes]]$table, annotation, by = "OG", all = T)
        }
    } else if (GeneDiscoveRobject$OrthologsTool == "OrthoFinder") {
        if (is.null(indexFilteredGenes)) {
            GeneDiscoveRobject$RunActive$N0Active <- merge(GeneDiscoveRobject$RunActive$N0Active, annotation, by.x = c("OG", "Gene Tree Parent Clade"), by.y = c("OG", "HOG"), all = T)
        } else {
            GeneDiscoveRobject$FilteredGenes[[indexFilteredGenes]]$table <- merge(GeneDiscoveRobject$FilteredGenes[[indexFilteredGenes]]$table, annotation, by.x = c("OG", "Gene Tree Parent Clade"), by.y = c("OG", "HOG"), all = T)
        }
    }

    return(GeneDiscoveRobject)
}

#' Maps gene annotations based on a GeneDiscoveR object and filtered genes index.
#'
#' This function takes a GeneDiscoveR object and an index of filtered genes as input,
#' and maps gene annotations based on the specified species with annotation.
#'
#' @param GeneDiscoveRobject A GeneDiscoveR object containing gene information.
#' @param indexFilteredGenes An index specifying the filtered genes to use for mapping annotations.
#' @param specieWithAnnotation The species with annotation to use for mapping.
#'
#' @return A tibble containing the mapped gene annotations.
#' @keywords internal
#' @examples
#' \dontrun{
#' # Map gene annotations for filtered genes at index 1
#' mapped_annotations <- .map_gene_annotation(GeneDiscoveRobject, indexFilteredGenes = 1)
#'
#' # Print the mapped annotations
#' print(mapped_annotations)
#' }
#' @importFrom readr read_tsv
#' @importFrom stringr str_split
#' @importFrom dplyr bind_rows rename
#' @export
.map_gene_annotation <- function(GeneDiscoveRobject = NULL, indexFilteredGenes = NULL, specieWithAnnotation = "MpTAKv6-Marchantia_polymorpha_rudelaris") {
    if (is.null(GeneDiscoveRobject$AnnotationFile) || !file.exists(GeneDiscoveRobject$AnnotationFile)) {
        stop("Error: 'GeneDiscoveRobject$AnnotationFile' is NULL or does not exist.")
    }
    tableAnnotation <- read_tsv(GeneDiscoveRobject$AnnotationFile, col_names = FALSE)
    Ids <- ""
    toAnnotation <- ""
    if (is.null(indexFilteredGenes)) {
        toAnnotation <- GeneDiscoveRobject$RunActive$N0Active
    } else {
        toAnnotation <- GeneDiscoveRobject$FilteredGenes[[indexFilteredGenes]]$table
    }
    result <- tibble()
    for (i in 1:dim(toAnnotation)[1]) {
        Ids <-
            str_split(toAnnotation[[specieWithAnnotation]][i], pattern = ", ")[[1]]
        OG <- toAnnotation$OG[i]
        if (GeneDiscoveRobject$OrthologsTool == "OrthoFinder") {
            GeneTree <- toAnnotation$`Gene Tree Parent Clade`[i]
        }
        count <- 0
        resultaux <- tibble()
        for (id in Ids) {
            if (!any(is.na(tableAnnotation[tableAnnotation[[1]] == str_split(id[1], pattern = "\\|")[[1]][1], 2]))) {
                if (GeneDiscoveRobject$OrthologsTool == "OrthoFinder") {
                    resultaux <- bind_rows(resultaux, tibble(
                        id, OG, GeneTree,
                        tableAnnotation[tableAnnotation[[1]] == str_split(id[1], pattern = "\\|")[[1]][1], -1]
                    ))
                } else if (GeneDiscoveRobject$OrthologsTool == "InParanoiDB") {
                    resultaux <- bind_rows(resultaux, tibble(
                        id, OG,
                        tableAnnotation[tableAnnotation[[1]] == str_split(id[1], pattern = "\\|")[[1]][1], -1]
                    ))
                }
                count <- count + 1
                if (count > 1) {
                    break
                }
            }
        }
        result <- bind_rows(result, resultaux)
    }
    if (GeneDiscoveRobject$OrthologsTool == "OrthoFinder") {
        result <- result %>% rename(
            Gene = id,
            OG = OG,
            HOG = GeneTree
        )
    } else if (GeneDiscoveRobject$OrthologsTool == "InParanoiDB") {
        result <- result %>% rename(
            Gene = id,
            OG = OG
        )
    }


    return(result)
}

#' Split annotation per gene
#'
#' This function takes a table of annotations for genes and splits the annotations
#' into separate columns for each gene. The annotations are assumed to be stored
#' as a character vector in the second column of the input table, with multiple
#' annotations separated by a specified separator.
#'
#' @param tableOG_HOG_Anno The input table of gene annotations.
#' @param sep The separator used to split the annotations (default is ";").
#'
#' @return A modified version of the input table with the annotations split into
#' separate columns for each gene.
#' @keywords internal
#' @examples
#' \dontrun{
#' # Create a sample table of gene annotations
#' tableOG_HOG_Anno <- data.frame(
#'     Gene = c("Gene1", "Gene2", "Gene3"),
#'     X2 = c(
#'         "Annotation1;Annotation2",
#'         "Annotation3",
#'         "Annotation4;Annotation5;Annotation6"
#'     )
#' )
#'
#' # Split the annotations per gene
#' result <- .split_annotation_per_gene(tableOG_HOG_Anno)
#' print(result)
#' }
#' @import dplyr
#' @export
.split_annotation_per_gene <- function(tableOG_HOG_Anno = NULL, sep = ";") {
    result <- tableOG_HOG_Anno %>%
        mutate(X2 = strsplit(as.character(X2), sep)) %>%
        unnest(X2) %>%
        group_by(Gene) %>%
        mutate(index = row_number()) %>%
        ungroup() %>%
        spread(index, X2, sep = "_")
    return(result)
}

#' Filter one gene per OG HOG
#'
#' This function filters a table of OG HOG annotations to keep only one gene per OG HOG combination.
#' The filtering is based on the number of annotations for each gene, keeping the gene with the minimum number of annotations.
#'
#' @param tableOG_HOG A data frame containing the OG HOG annotations.
#'
#' @return A data frame with one gene per OG HOG combination, based on the minimum number of annotations.
#' @keywords internal
#' @examples
#' \dontrun{
#' # Create a sample table of OG HOG annotations
#' tableOG_HOG <- data.frame(
#'     OG = c("OG1", "OG1", "OG2", "OG2", "OG3"),
#'     HOG = c("HOG1", "HOG1", "HOG2", "HOG2", "HOG3"),
#'     Gene = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5"),
#'     Annotation1 = c(1, 1, 1, 1, 1),
#'     Annotation2 = c(1, 1, 1, 1, 1),
#'     Annotation3 = c(1, 1, 1, 1, 1)
#' )
#'
#' # Filter one gene per OG HOG
#' filtered_table <- .filter_oneGene_per_OG_HOG(tableOG_HOG)
#' }
#' @import dplyr
#' @export
.filter_oneGene_per_OG_HOG <- function(tableOG_HOG = NULL) {
    result <- tableOG_HOG %>%
        mutate(Nannotations = rowSums(!is.na(.[-(1:3)]))) %>%
        group_by(OG, HOG, Gene) %>%
        summarise(Nannotations = min(Nannotations), .groups = "drop") %>%
        ungroup() %>%
        group_by(OG, HOG) %>%
        arrange(Nannotations) %>%
        slice(1) %>%
        ungroup() %>%
        left_join(tableOG_HOG, by = c("OG", "HOG", "Gene"))
    return(result)
}
