#' Merge inParanoid Files
#'
#' This function merges multiple inParanoid files into a single dataframe.
#'
#' @param directory The directory where the inParanoid files are located.
#'
#' @return A merged dataframe containing the data from all the inParanoid files.
#' @keywords internal
#' @examples
#' # Set the directory where the inParanoid files are located
#' directory <- system.file("extdata", "inparanoidbFiles", package = "GeneDiscoveR")
#'
#' # Merge the inParanoid files
#' merged_df <- .merge_inParanoid_files(directory)
#'
#' # Print the merged dataframe
#' \dontrun{
#' print(merged_df)
#' }
#'
#' @import readr dplyr stringr
#' @export
.merge_inParanoid_files <- function(directory) {
    options(readr.show_progress = FALSE)
    file_list <- list.files(directory, pattern = "^inParanoid.*", full.names = TRUE)

    dfs <- lapply(file_list, function(file) {
        df <- suppressWarnings(suppressMessages(read_tsv(file)))
        df <- df %>% mutate(OrtoBsplit = str_split(OrtoB, "\t"))
        return(df)
    })

    dfMerge <- Reduce(function(x, y) full_join(x, y, by = "OrtoA"), dfs)

    return(dfMerge)
}

#' Check Orthologs in Ortholog Set B
#'
#' This function checks for orthologs in Ortholog Set B based on the values in Ortholog Set A.
#'
#' @param dfMerge A data frame containing the merged data from Ortholog Set A and Ortholog Set B.
#' @param cores The number of cores to use for parallel processing.
#' @return A list of unique indices of orthologs in Ortholog Set B that match the values in Ortholog Set A.
#' @keywords internal
#' @examples
#' df <- data.frame(OrtoA = c("A", "B", "C"), OrtoB.1 = c("A", "B", "D"), OrtoB.2 = c("B", "C", "E"))
#' .check_OrtoA_in_OrtoB(df, cores = 4)
#' # Output: [[1]]
#' # [1] 1 2
#' #
#' # [[2]]
#' # [1] 2 3
#' #
#' # [[3]]
#' # integer(0)
#' @import foreach
#' @import doParallel
#' @import stringr
#' @import parallel
#' @keywords internal
#' @export
.check_OrtoA_in_OrtoB <- function(dfMerge, cores = 1) {
    cl <- makeCluster(cores)
    registerDoParallel(cl)

    ortoB_cols <- grep("^OrtoB\\.", names(dfMerge), value = TRUE)

    resultados <- foreach(i = seq_along(dfMerge$OrtoA), .combine = "c", .packages = "stringr") %dopar% {
        ortoA_val <- dfMerge$OrtoA[i]
        coincidencias <- integer(0)

        for (col in ortoB_cols) {
            coincidencias <- c(coincidencias, which(str_detect(dfMerge[[col]], ortoA_val)))
        }

        list(unique(coincidencias))
    }
    stopCluster(cl)

    return(resultados)
}

#' Get the even values from a vector
#'
#' This function takes a vector as input and returns a new vector containing only the even values from the input vector.
#'
#' @param x A numeric vector
#' @return A numeric vector containing only the even values from the input vector
#' @keywords internal
#' @examples
#' x <- c(1, 2, 3, 4, 5, 6)
#' .get_even_values(x)
#'
#' # Output: 2 4 6
#' @export
.get_even_values <- function(x) {
    odd_indices <- seq_along(x) %% 2 == 0
    x[odd_indices]
}

#' Merge groups based on given results
#'
#' This function merges groups based on the given results. It takes a data frame
#' \code{dfMerge} and a list \code{resultados} as input. It returns a list
#' containing two data frames: \code{dfMergemodifyunprocess} and
#' \code{dfMergemodifyprocess}.
#'
#' @param dfMerge A data frame containing the data to be merged
#' @param resultados A list of results used for merging groups
#' @keywords internal
#' @return A list containing two data frames: \code{dfMergemodifyunprocess} and
#' \code{dfMergemodifyprocess}
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'     OrtoA = c("A", "B", "C"),
#'     OrtoBsplit1 = c("1,2,3", "4,5,6", "7,8,9"),
#'     OrtoBsplit2 = c("10,11,12", "13,14,15", "16,17,18")
#' )
#' resultados <- list(c("1\t2\t3", "4\t5\t6"), c("7\t8\t9"))
#' .merge_groups(df, resultados)
#' }
#' @export
.merge_groups <- function(dfMerge, resultados) {
    dfMergemodifyunprocess <- c()
    dfMergemodifyprocess <- c()
    ortoB_cols <- grep("^OrtoBsplit", names(dfMerge), value = TRUE)
    indices <- c()
    success <- c()
    for (i in seq_along(dfMerge$OrtoA)) {
        if (length(resultados[[i]]) > 0 && !any(resultados[[i]] %in% success)) {
            indices <- c(indices, i)
            for (k in 1:length(unlist(str_split(resultados[[i]], "\t")))) {
                for (j in seq(i + 1, length(resultados))) {
                    if (unlist(str_split(resultados[[i]], "\t"))[k] %in% paste0(unlist(str_split(resultados[[j]], "\t")), collapse = "")) {
                        indices <- c(indices, j)
                    }
                }
            }
            indices <- unique(c(indices, unlist(str_split(resultados[[i]], "\t"))))
            success <- unique(c(success, indices, unlist(str_split(resultados[[i]], "\t"))))

            ids <- vector("list", length(indices))
            for (j in seq_along(indices)) {
                for (col in ortoB_cols) {
                    ids[[j]] <- c(ids[[j]], map(dfMerge[indices[j], col][[1]], ~ .get_even_values(eval(.x))))
                }
            }

            lista_combinada <- vector("list", length(ortoB_cols))
            for (j in 1:length(ortoB_cols)) {
                elementos_a_unir <- lapply(ids, function(x) if (length(x) >= j) x[[j]] else numeric(0))
                if (!is.null(unique(unlist(elementos_a_unir)))) {
                    lista_combinada[[j]] <- unique(unlist(elementos_a_unir))
                } else {
                    lista_combinada[[j]] <- numeric(0)
                }
            }
            renglon <- sapply(lista_combinada, function(x) paste0(unlist(x), collapse = ","))
            renglon <- as.data.frame(t(renglon), stringsAsFactors = FALSE)
            colnames(renglon) <- ortoB_cols
            dfMergemodifyprocess <- rbind(dfMergemodifyprocess, cbind(dfMerge$OrtoA[i], renglon))
        } else if (!i %in% success) {
            dfMergemodifyunprocess <- rbind(dfMergemodifyunprocess, cbind(dfMerge$OrtoA[i], dfMerge[i, ortoB_cols]))
        }
        indices <- c()
    }
    return(list(dfMergemodifyunprocess = dfMergemodifyunprocess, dfMergemodifyprocess = dfMergemodifyprocess))
}

#' Merge groups in parallel
#'
#' This function merges groups in parallel using the provided dataframe.
#'
#' @param dfMergemodify The dataframe to be modified and merged.
#' @param cores The number of cores to be used for parallel processing. Default is 1.
#' @return The modified and merged dataframe.
#' @keywords internal
#' @import dplyr
#' @import foreach
#' @import doParallel
#' @import parallel
#' @importFrom purrr map
#' @export
.merge_groups_parallel <- function(dfMergemodify, cores = 1) {
    ortoB_cols <- grep("^OrtoBsplit", names(dfMergemodify), value = TRUE)

    cl <- makeCluster(cores)

    registerDoParallel(cl)

    resultados1 <- foreach(i = seq_along(dfMergemodify$`dfMerge$OrtoA[i]`), .packages = c("dplyr", "purrr"), .export = ".get_even_values") %dopar% {
        fila <- dfMergemodify[i, ]
        for (col in ortoB_cols) {
            fila <- fila %>% mutate(!!col := map(fila[col][[1]], ~ .get_even_values(eval(.x))))
        }
        fila
    }
    dfMerge_paralelo <- bind_rows(resultados1)
    stopCluster(cl)
    dfMergemodify <- dfMerge_paralelo

    return(dfMergemodify)
}

#' Update the dfMerge dataframe with orthologous gene IDs from InParanoid database
#'
#' This function takes a dataframe (dfMerge) and updates it by replacing the orthologous gene IDs
#' in the columns starting with "OrtoBsplit" with cleaned IDs obtained from the InParanoid database.
#' The function requires the following parameters:
#'
#' \code{dfMerge}: The input dataframe to be updated.
#'
#' \code{at_prefix}: The prefix used in the gene IDs to be replaced.
#'
#' \code{datafile}: The datafile containing the orthologous gene IDs from the InParanoid database.
#'
#' \code{columnID}: The column name in the datafile that contains the orthologous gene IDs.
#'
#' The function iterates over each row of the dfMerge dataframe and for each row, it extracts the
#' orthologous gene IDs from the columns starting with "OrtoBsplit". It then cleans the IDs by removing
#' the prefix and any additional characters. The cleaned IDs are then replaced in the dfMerge dataframe.
#' Additionally, the function performs some data cleaning operations on the dfMerge dataframe, such as
#' removing leading and trailing whitespaces and removing consecutive commas in the columns starting with "OrtoBsplit".
#' Finally, the function renames the columns in the dfMerge dataframe based on the orthologous gene IDs
#' obtained from the datafile.
#'
#' @param dfMerge The input dataframe to be updated.
#' @param at_prefix The prefix used in the gene IDs to be replaced.
#' @param datafile The datafile containing the orthologous gene IDs from the InParanoid database.
#' @param columnID The column name in the datafile that contains the orthologous gene IDs.
#'
#' @return The updated dfMerge dataframe with cleaned orthologous gene IDs.
#' @keywords internal
#' @examples
#' # Example usage of the update_dfMerge function
#' \dontrun{
#' df <- .update_dfMerge(dfMerge, "AT", datafile, "geneID")
#' }
#' @export
.update_dfMerge <- function(dfMerge, at_prefix, datafile, columnID) {
    dfMerge1 <- dfMerge
    ortoB_cols <- grep("^OrtoBsplit", names(dfMerge), value = TRUE)
    ids_ortoB_limpio <- c()
    for (i in seq_along(dfMerge1$OrtoA)) {
        ortoA_actual <- dfMerge1$OrtoA[i]

        for (col in ortoB_cols) {
            ids_ortoB <- eval(dfMerge1[i, col])[[1]]
            ids_AT <- str_extract_all(ids_ortoB, paste0(at_prefix, "[^,]+")) %>% unlist()
            ids_ortoB_limpio <- str_replace_all(ids_ortoB, paste0(at_prefix, "[^,]+"), "") %>% unlist()
            ids_ortoB_limpio <- eval(ids_ortoB_limpio)
            ids_ortoB_limpio <- ids_ortoB_limpio[nzchar(ids_ortoB_limpio)]
            ids_ortoB_limpio <- paste(ids_ortoB_limpio, collapse = ", ")

            if (length(ids_AT) > 0) {
                check <- str_detect(pattern = unique(ids_AT), string = ortoA_actual)
            }

            if (!length(ids_AT) == 0 && !all(check)) {
                ortoA_actual <- paste(ortoA_actual, unique(ids_AT[!check]), sep = ", ", collapse = ", ")
            }

            dfMerge1[i, col] <- ids_ortoB_limpio
        }

        dfMerge1$OrtoA[i] <- ortoA_actual
    }
    dfMerge1 <- dfMerge1 %>%
        mutate(across(starts_with("OrtoBsplit"), ~ str_trim(.x))) %>%
        mutate(across(starts_with("OrtoBsplit"), ~ str_replace_all(.x, ",{2,}", ""))) %>%
        mutate(across(starts_with("OrtoBsplit"), ~ str_replace_all(.x, "^(,)+", "")))
    inparanoidbID <- sort(datafile[[columnID]])
    dfMerge1 <- dfMerge1 %>%
        rename_with(~ inparanoidbID[-which(principalSpecie %in% inparanoidbID)], starts_with("OrtoBsplit")) %>%
        rename(!!sym(inparanoidbID[which(principalSpecie %in% inparanoidbID)]) := OrtoA)

    return(dfMerge1)
}
