set_annotation_file <- function(PhenoRobject = NULL, annotationFile = NULL) {
    if (is.null(PhenoRobject) || is.null(annotationFile)) {
        stop("Error: 'PhenoRobject' and 'annotationFile' cannot be NULL.")
    }
    PhenoRobject$AnnotationFile <- annotationFile
    return(PhenoRobject)
}
map_annotation <- function(PhenoRobject = NULL, indexFilteredGenes = NULL, specieWithAnnotation = "MpTAKv6-Marchantia_polymorpha_rudelaris") {
    annotation <- map_gene_annotation(PhenoRobject = PhenoRobject, indexFilteredGenes = indexFilteredGenes, specieWithAnnotation = specieWithAnnotation)
    annotation <- split_annotation_per_gene(tableOG_HOG_Anno = annotation)
    annotation <- filter_oneGene_per_OG_HOG(tableOG_HOG = annotation)
    PhenoRobject$FilteredGenes[[indexFilteredGenes]]$table <- merge(PhenoRobject$FilteredGenes[[indexFilteredGenes]]$table, annotation, by.x = c("OG", "Gene Tree Parent Clade"), by.y = c("OG", "HOG"), all = T)
    return(PhenoRobject)
}

map_gene_annotation <- function(PhenoRobject = NULL, indexFilteredGenes = NULL, specieWithAnnotation = "MpTAKv6-Marchantia_polymorpha_rudelaris") {
    if (is.null(PhenoRobject$AnnotationFile) || !file.exists(PhenoRobject$AnnotationFile)) {
        stop("Error: 'PhenoRobject$AnnotationFile' is NULL or does not exist.")
    }
    tableAnnotation <- read_tsv(PhenoRobject$AnnotationFile, col_names = FALSE)
    # Check if the tableAnnotation is read correctly
    result <- c()
    Ids <- ""
    for (i in 1:dim(PhenoRobject$FilteredGenes[[indexFilteredGenes]]$table)[1]) {
        Ids <-
            str_split(PhenoRobject$FilteredGenes[[indexFilteredGenes]]$table[[specieWithAnnotation]][i], pattern = ", ")[[1]]
        OG <- PhenoRobject$FilteredGenes[[indexFilteredGenes]]$table$OG[i]
        GeneTree <- PhenoRobject$FilteredGenes[[indexFilteredGenes]]$table$`Gene Tree Parent Clade`[i]
        for (id in Ids) {
            if (!any(is.na(tableAnnotation[tableAnnotation$X1 == str_split(id[1], pattern = "\\|")[[1]][1], 2]))) {
                result <- rbind(result, c(
                    id, OG, GeneTree,
                    tableAnnotation[tableAnnotation$X1 == str_split(id[1], pattern = "\\|")[[1]][1], 2]
                ))
            } else {
                result <- rbind(result, c("NAgene", OG, GeneTree, "NAannotation"))
            }
        }
    }
    result <- as.data.frame(result)
    result <- result %>% rename(
        Gene = V1,
        OG = V2,
        HOG = V3,
        Anno = X2
    )
    return(result)
}
split_annotation_per_gene <- function(tableOG_HOG_Anno = NULL) {
    result <- tableOG_HOG_Anno %>%
        mutate(Anno = strsplit(as.character(Anno), ";")) %>%
        unnest(Anno) %>%
        group_by(Gene) %>%
        mutate(index = row_number()) %>%
        ungroup() %>%
        spread(index, Anno, sep = "_")
    return(result)
}
filter_oneGene_per_OG_HOG <- function(tableOG_HOG = NULL) {
    result <- tableOG_HOG %>%
        mutate(Nannotations = rowSums(!is.na(.[-(1:2)]))) %>%
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
