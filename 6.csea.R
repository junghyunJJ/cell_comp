rm(list = ls())
library(tidyverse)
library(data.table)
library(tidytext)

library(jjutil)
library(openxlsx)

library(fgsea)

rawres <- readRDS("res/summary_entropy.rds")
res <- rawres %>% 
    mutate(cellgroup = if_else(grepl(";", annotation), "hetero", "homo"))

sel <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# sel <- 1
types <- grep("entropy|genescore", colnames(res), value = TRUE)

types <- types[c(1:4, 9, 10, 25, 26, 29, 30, 35, 36)]
#  [1] "win_entropy"                        "win_genescore"                     
#  [3] "lose_entropy"                       "lose_genescore"                    
#  [5] "apoptosis_msigdbhallmark_entropy"   "apoptosis_msigdbhallmark_genescore"
#  [7] "HALLMARK_MYC_TARGETS_V1_entropy"    "HALLMARK_MYC_TARGETS_V1_genescore" 
#  [9] "HALLMARK_P53_PATHWAY_entropy"       "HALLMARK_P53_PATHWAY_genescore"    
# [11] "winlose_diff_entropy"               "winlose_diff_genescore"     
sel_type <- types[sel]

#########################################################################################
#########################################################################################
#########################################################################################
list_rank <- list()
list_cellset <- list()
list_fgsea <- list()

set.seed(1)

# for (sel_time in c("E9.5", "E10.5")) {
for (sel_time in c("E9.5", "E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16.5")) {
    # sel_time <- "E9.5"
    print(sel_time)
    dat <- res %>%
        filter(timepoint == sel_time) %>%
        filter(!is.na(get(sel_type)))

    # rank data
    rawrankdat <- dat %>%
        select(cellname, get("sel_type")) %>%
        arrange(get(sel_type))

    rankdat <- rawrankdat[, 2]
    names(rankdat) <- rawrankdat$cellname

    # hetero geneset data
    hetero_dat <- dat %>%
        filter(grepl(";", annotation)) %>%
        mutate(cellset_id = str_glue("{sel_time}_{annotation}_{celltype}_{n_grid}"))

    tot_hetero_cellset <- hetero_dat$cellset_id %>% unique()
    cellset <- list()
    
    for (sel_hetero_cellset in tot_hetero_cellset) {
        # sel_hetero_cellset <- tot_hetero_cellset[1]
        print(sel_hetero_cellset)
        sel_hetero_cells <- hetero_dat %>%
            filter(cellset_id == sel_hetero_cellset) %>%
            select(cellname) %>%
            pull()

        cellset <- append(cellset, list(sel_hetero_cells))
    }
    names(cellset) <- tot_hetero_cellset
    
    # run CSEA
    fgseaRes <- fgseaMultilevel(pathways = cellset, stats = rankdat)

    # save
    list_fgsea <- append(list_fgsea, list(fgseaRes))
    list_rank <- append(list_rank, list(rankdat))
    list_cellset <- append(list_cellset, cellset)
}


names(list_rank) <- paste0(sel_type, "-", c("E9.5", "E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16.5"))
# names(list_fgsea) <- names(list_rank) <- paste0(sel_type, "-", c("E9.5", "E10.5"))

save(list_fgsea, list_rank, list_cellset, file = paste0("res/res_csea/", sel_type, ".RData"))

# list_fgsea[[1]] %>% arrange(padj) %>% filter(padj < 0.05, ES > 0)
# list_fgsea[[2]] %>% arrange(padj) %>% filter(padj < 0.05, ES > 0)
# names(list_fgsea)
