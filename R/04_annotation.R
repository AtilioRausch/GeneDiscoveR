set_annotation_file <- function(GeneDiscoveRobject = NULL, annotationFile = NULL) {
    if (is.null(GeneDiscoveRobject) || is.null(annotationFile)) {
        stop("Error: 'GeneDiscoveRobject' and 'annotationFile' cannot be NULL.")
    }
    GeneDiscoveRobject$AnnotationFile <- annotationFile
    return(GeneDiscoveRobject)
}
map_annotation <- function(GeneDiscoveRobject = NULL, indexFilteredGenes = NULL, oneColumn = TRUE, specieWithAnnotation = "MpTAKv6-Marchantia_polymorpha_rudelaris") {
    annotation <- map_gene_annotation(GeneDiscoveRobject = GeneDiscoveRobject, indexFilteredGenes = indexFilteredGenes, specieWithAnnotation = specieWithAnnotation)
    if (oneColumn) {
        annotation <- split_annotation_per_gene(tableOG_HOG_Anno = annotation)
    }
    annotation <- filter_oneGene_per_OG_HOG(tableOG_HOG = annotation)
    GeneDiscoveRobject$FilteredGenes[[indexFilteredGenes]]$table <- merge(GeneDiscoveRobject$FilteredGenes[[indexFilteredGenes]]$table, annotation, by.x = c("OG", "Gene Tree Parent Clade"), by.y = c("OG", "HOG"), all = T)
    return(GeneDiscoveRobject)
}
map_gene_annotation <- function(GeneDiscoveRobject = NULL, indexFilteredGenes = NULL, specieWithAnnotation = "MpTAKv6-Marchantia_polymorpha_rudelaris") {
    if (is.null(GeneDiscoveRobject$AnnotationFile) || !file.exists(GeneDiscoveRobject$AnnotationFile)) {
        stop("Error: 'GeneDiscoveRobject$AnnotationFile' is NULL or does not exist.")
    }
    tableAnnotation <- read_tsv(GeneDiscoveRobject$AnnotationFile, col_names = FALSE)
    # Check if the tableAnnotation is read correctly
    Ids <- ""
    result <- tibble()
    for (i in 1:dim(GeneDiscoveRobject$FilteredGenes[[indexFilteredGenes]]$table)[1]) {
        Ids <-
            str_split(GeneDiscoveRobject$FilteredGenes[[indexFilteredGenes]]$table[[specieWithAnnotation]][i], pattern = ", ")[[1]]
        OG <- GeneDiscoveRobject$FilteredGenes[[indexFilteredGenes]]$table$OG[i]
        GeneTree <- GeneDiscoveRobject$FilteredGenes[[indexFilteredGenes]]$table$`Gene Tree Parent Clade`[i]
        count <- 0
        resultaux <- tibble()
        for (id in Ids) {
            if (!any(is.na(tableAnnotation[tableAnnotation[[1]] == str_split(id[1], pattern = "\\|")[[1]][1], 2]))) {
                resultaux <- bind_rows(resultaux, tibble(
                    id, OG, GeneTree,
                    tableAnnotation[tableAnnotation[[1]] == str_split(id[1], pattern = "\\|")[[1]][1], -1]
                ))
                count <- count + 1
                if (count > 1) {
                    break
                }
            }
        }
        result <- bind_rows(result, resultaux)
    }
    result <- result %>% rename(
        Gene = id,
        OG = OG,
        HOG = GeneTree
    )
    return(result)
}
split_annotation_per_gene <- function(tableOG_HOG_Anno = NULL) {
    # ";"
    sep <- "\t"
    result <- tableOG_HOG_Anno %>%
        mutate(Anno = strsplit(as.character(Anno), sep)) %>%
        unnest(Anno) %>%
        group_by(Gene) %>%
        mutate(index = row_number()) %>%
        ungroup() %>%
        spread(index, Anno, sep = "_")
    return(result)
}
filter_oneGene_per_OG_HOG <- function(tableOG_HOG = NULL) {
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
