rm(list = ls())
library(tidyverse)
library(data.table)

library(jjutil)
library(tidytext)

scaling <- 1
rawmosta <- readRDS("res/new_summary_mosta.rds")
mosta <- rawmosta %>%
    mutate(winlose_diff = (win - lose) / scaling) %>%
    mutate(proapop_diff = (proliferation_plos - apoptosis_msigdbhallmark) / scaling) %>%
    mutate(winapop_diff = (win - apoptosis_msigdbhallmark) / scaling) %>% 
    mutate(winloseapop_diff = (win - apoptosis_msigdbhallmark - lose) / scaling)

# value_names <- colnames(mosta)[c(6, 12:ncol(mosta))]
# t_mosta <- cbind(
#     mosta[, !colnames(mosta) %in% value_names],
#     exp(mosta[, colnames(mosta) %in% value_names])
# )

#################################################################################
### Cor of Marker scores ########################################################
#################################################################################

# seltime <- "E9.5"
# seltissue <- "Brain"

lapply(levels(mosta$timepoint), function(seltime) {  
    sel_mosta <- mosta %>%
        filter(timepoint == seltime)
    
    # all tisses
    sel_mosta %>%
        select(-cell_name, -newy, -newx, -x, -y, -breaks_ratio, -annotation, -timepoint, -section, -time) %>% 
        cor %>% 
        # ggcorrplot::ggcorrplot(hc.order = T, type = "lower", lab = TRUE) + 
        ggcorrplot::ggcorrplot(lab = TRUE) +
        ggtitle(str_glue("{seltime}"))
    ggsave(str_glue("fig/cor/cor_{seltime}.png"), width = 30, height = 40, units = "in", bg = "white")

    # separate tissue
    pbmcapply::pbmclapply(unique(sel_mosta$annotation), function(seltissue) {
        #cat(str_glue("{seltime}_{seltissue}"), "\n")
        sel_mosta %>%
            filter(annotation == seltissue) %>% 
            select(-cell_name, -newy, -newx, -x, -y, -breaks_ratio, -annotation, -timepoint, -section, -time) %>% 
            cor %>% 
            # ggcorrplot::ggcorrplot(hc.order = T, type = "lower", lab = TRUE) + 
            ggcorrplot::ggcorrplot(lab = TRUE) +
            ggtitle(str_glue("{seltime}_{seltissue}"))
        ggsave(str_glue("fig/cor/cor_{seltime}_{seltissue}.png"), width = 30, height = 40, units = "in", bg = "white")
    }, mc.cores = 10)
})


#################################################################################
### Shannon entropy of cell marker score ########################################
#################################################################################

selfiles <- list.files("res/res_entropy/", full.names = TRUE, pattern = str_glue("tissue"))
res_tissue <- pbmcapply::pbmclapply(selfiles, function(sel) {
    readRDS(sel)
}, mc.cores = 20) %>% rbindlist()
res_tissue$analysis <- "tissue"
saveRDS(res_tissue, "res/summary_entropy_tissue.rds")


lapply(2:5, function(cr) {
    cat(str_glue("window_{cr}"), "\n")
    selfiles <- list.files("res/res_entropy/", full.names = TRUE, pattern = str_glue("window_{cr}"))

    res <- pbmcapply::pbmclapply(selfiles, function(sel) {
        readRDS(sel)
    }, mc.cores = 20) %>% rbindlist()
    res$analysis <- str_glue("window_{cr}")
    saveRDS(res, str_glue("res/summary_entropy_window_{cr}.rds"))
})





# rawres <- readRDS("res/summary_entropy_FALSE_2.rds")
# res <- rawres %>% 
#     filter(type %in% c("ratio", "winlose_diff", "proapop_diff", "winapop_diff", "winloseapop_diff",
#         "win", "lose", "proliferation_plos", "apoptosis_msigdbhallmark", "essential_both", "haploinsufficiency",
#         "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"
#     )) %>% 
#     mutate(timepoint = factor(timepoint, levels = c("E9.5", "E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16.5"))) %>% 
#     mutate(col = if_else(annotation %in% c("Brain", "Spinal cord"), "1", "0"))
# res$type %>% factor %>% levels


# res$type %>% factor %>% levels
# factor(res$type, levels = c(
#     "ratio", "winlose_diff", "proapop_diff", "winapop_diff", "winloseapop_diff",
#     "win", "lose", "proliferation_plos", "apoptosis_msigdbhallmark", "essential_both", "haploinsufficiency",
#     "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"
# )) %>% levels()
# res$type <- factor(res$type, levels = c(
#         "ratio", "winlose_diff", "proapop_diff", "winapop_diff", "winloseapop_diff",
#         "win", "lose", "proliferation_plos", "apoptosis_msigdbhallmark", "essential_both", "haploinsufficiency",
#         "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"
#     ), 
#     labels = c(
#         "ratio", "winlose_diff", "proapop_diff", "winapop_diff", "winloseapop_diff",
#         "win", "lose", "proliferation", "apoptosis", "essential", "haploinsufficiency",
#         "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"))

# sel <- "ratio"
# lapply(levels(res$type), function(sel) {
#     cat(sel, "\n")
#     # sel <- "win"
#     if (sel %in% c("ratio", "winlose_diff", "proapop_diff", "winapop_diff", "winloseapop_diff", "win", "lose")) {
#         set_col <- c("white", "firebrick3")
#     } else {
#         set_col <- c("white", "royalblue3")
#     }
#     res %>% 
#         filter(type == sel) %>%
#         ggplot(aes(entropy, reorder_within(annotation, entropy, timepoint, median), fill = col)) +
#         theme_bw() +
#         geom_boxplot() +
#         scale_fill_manual(values = set_col) +
#         facet_wrap(~ timepoint, scales = "free", nrow = 4) +
#         theme(legend.position = "none") +
#         labs(y = "", title = sel)
# r    ggsave(str_glue("fig/entropy/entropy_{sel}.png"), width = 15, height = 20, units = "in", bg = "white")
# })