rm(list = ls())
library(tidyverse)
library(data.table)
library(tidytext)

library(jjutil)
library(openxlsx)

# library(Vennerable)
library(ggvenn)

library(fgsea)

#  [1] "win_entropy"                        "win_genescore"
#  [3] "lose_entropy"                       "lose_genescore"
#  [5] "apoptosis_msigdbhallmark_entropy"   "apoptosis_msigdbhallmark_genescore"
#  [7] "HALLMARK_MYC_TARGETS_V1_entropy"    "HALLMARK_MYC_TARGETS_V1_genescore"
#  [9] "HALLMARK_P53_PATHWAY_entropy"       "HALLMARK_P53_PATHWAY_genescore"
# [11] "winlose_diff_entropy"               "winlose_diff_genescore"

timepoint <- c("E9.5", "E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16.5")

list_res_cellgroup <- list_res_cell <- list_gg_res_cellgroup <- list_gg_res_cell <- list()


for (i in seq_len(length(timepoint))) {
    print(timepoint[i])
    # i <- 4
    # sel_type <- "win_genescore"
    # load(paste0("res/res_csea/", sel_type, ".RData"))
    # win <- list_fgsea[[i]] %>% arrange(pval) %>% filter(padj < 0.05, NES < 0)

    # sel_type <- "lose_genescore"
    # load(paste0("res/res_csea/", sel_type, ".RData"))
    # lose <- list_fgsea[[i]] %>% arrange(pval) %>% filter(padj < 0.05, NES > 0)

    sel_type <- "apoptosis_msigdbhallmark_genescore"
    load(paste0("res/res_csea/", sel_type, ".RData"))
    apop <- list_fgsea[[i]] %>% arrange(pval) %>% filter(padj < 0.05, NES > 0)

    sel_type <- "HALLMARK_P53_PATHWAY_genescore"
    load(paste0("res/res_csea/", sel_type, ".RData"))
    p53 <- list_fgsea[[i]] %>% arrange(pval) %>% filter(padj < 0.05, NES > 0)
    
    sel_type <- "winlose_diff_genescore"
    load(paste0("res/res_csea/", sel_type, ".RData"))
    comp_rate <- list_fgsea[[i]] %>% arrange(pval) %>% filter(padj < 0.05, NES < 0)

    sel_type <- "HALLMARK_MYC_TARGETS_V1_genescore"
    load(paste0("res/res_csea/", sel_type, ".RData"))
    myc <- list_fgsea[[i]] %>% arrange(pval) %>% filter(padj < 0.05, NES < 0)

    x <- list(
        apoptosis = apop$pathway,
        `P53 pathway\n(HALLMARK)` = p53$pathway,
        `MYC targets V1\n(HALLMARK)` = myc$pathway,
        `competition rate` = comp_rate$pathway
    )
    gg <- ggvenn::ggvenn(x, fill_color = c("orange", "pink", "white", "yellow"), stroke_size = 0.2, set_name_size = 3, show_outside = "always") +
        labs(title = str_glue("{timepoint[i]}: cell group (hetero data)"))
    list_gg_res_cellgroup <- append(list_gg_res_cellgroup, list(gg))

    v_table <- Vennerable::Venn(x)
    sel_cellgroup <- c(v_table@IntersectionSets$`1111`)

    # leadingEdge cells
    sel_apop <- apop %>%
        filter(pathway %in% sel_cellgroup) %>%
        select(leadingEdge) %>% pull %>% unlist

    sel_p53 <- p53 %>%
        filter(pathway %in% sel_cellgroup) %>%
        select(leadingEdge) %>%
        pull() %>%
        unlist()

    sel_comp_rate <- comp_rate %>%
        filter(pathway %in% sel_cellgroup) %>%
        select(leadingEdge) %>%
        pull() %>%
        unlist()
    
    sel_myc <- myc %>%
        filter(pathway %in% sel_cellgroup) %>%
        select(leadingEdge) %>%
        pull() %>%
        unlist()
    
    xx <- list(
        apoptosis = sel_apop,
        `P53 pathway\n(HALLMARK)` = sel_p53,
        `MYC targets V1\n(HALLMARK)` = sel_comp_rate,
        `competition rate` = sel_myc
    )
    ggvenn::ggvenn(xx)
    sel_v_table <- Vennerable::Venn(xx)
    sel_cell <- c(sel_v_table@IntersectionSets$`1111`)

    gg2 <- ggvenn::ggvenn(xx, fill_color = c("orange", "pink", "white", "yellow"), stroke_size = 0.2, set_name_size = 3, show_outside = "always") +
        labs(title = str_glue("{timepoint[i]}: cell (hetero data)"))
    list_gg_res_cell <- append(list_gg_res_cell, list(gg2))

    list_res_cellgroup <- append(list_res_cellgroup, list(sel_cellgroup))
    list_res_cell <- append(list_res_cell, list(sel_cell))
}

names(list_res_cellgroup) <- names(list_res_cell) <- names(list_gg_res_cellgroup) <- names(list_gg_res_cell) <- timepoint
save(list_res_cellgroup, list_res_cell, list_gg_res_cellgroup, list_gg_res_cell, file = "res/csea_genescore.Rdata")
