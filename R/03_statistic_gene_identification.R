gene_identification_by_phenotype <- function(formula, PhenoRobject = NULL, statistic = "Fisher", name = "PerType") {
    if (statistic == "Fisher") {
        # Comprueba si ambos elementos de la fórmula están en PhenoRobject$Phenotypes
        formula_elements <- c(formula[[2]], formula[[3]])
        if (!(all(formula_elements %in% names(PhenoRobject$Phenotypes)))) {
            stop("Both elements of the formula must be elements of PhenoRobject$Phenotypes.")
        } else {
            nameColumn1 <- paste0("N_", str_replace_all(pattern = " ", replacement = "_", formula[[2]]))
            nameColumn2 <- paste0("N_", str_replace_all(pattern = " ", replacement = "_", formula[[3]]))
            # Realiza la prueba exacta de Fisher
            PhenoRobject$RunActive$N0Active <- PhenoRobject$RunActive$N0Active %>%
                rowwise() %>%
                mutate(
                    nameColumn1 = sum(!is.na(c_across(all_of(PhenoRobject$Phenotypes[[formula[[2]]]])))),
                    nameColumn2 = sum(!is.na(c_across(all_of(PhenoRobject$Phenotypes[[formula[[3]]]]))))
                )
            nameColumn3 <- paste0("fisherResult", name)
            PhenoRobject$RunActive$N0Active <- PhenoRobject$RunActive$N0Active %>%
                rowwise() %>%
                mutate(
                    !!nameColumn3 := list(fisher.test(matrix(
                        c(
                            (nameColumn1),
                            (nameColumn2),
                            length(PhenoRobject$Phenotypes[[formula[[2]]]]) - (nameColumn1),
                            length(PhenoRobject$Phenotypes[[formula[[3]]]]) - (nameColumn2)
                        ),
                        nrow = 2, byrow = TRUE
                    )))
                ) %>%
                rowwise() %>%
                mutate(
                    !!paste0("pvalueFisher", name) := (!!sym(nameColumn3))[["p.value"]],
                    !!paste0("oddsRatioFisher", name) := (!!sym(nameColumn3))[["estimate"]]
                ) %>%
                rename(!!nameColumn1 := nameColumn1, !!nameColumn2 := nameColumn2)
        }
        return(PhenoRobject)
    }
}
