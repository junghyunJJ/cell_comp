rm(list = ls())
library(tidyverse)
library(data.table)

library(jjutil)
library(tidytext)

#################################################################################
### summary of entropy resutls ##################################################
#################################################################################

# mosta <- readRDS("res/summary_mosta.rds")

# # set param file for Slurm array
# tmp_type <- grep("_entropy$", colnames(mosta), value = TRUE)
# (type <- sub("_entropy", "", tmp_type))
# #  [1] "win"                                 "lose"                               
# #  [3] "proliferation_moscot"                "proliferation_plos"                 
# #  [5] "apoptosis_msigdbhallmark"            "essential_all"                      
# #  [7] "essential_both"                      "haploinsufficiency"                 
# #  [9] "HALLMARK_WNT_BETA_CATENIN_SIGNALING" "HALLMARK_IL6_JAK_STAT3_SIGNALING"   
# # [11] "HALLMARK_MTORC1_SIGNALING"           "HALLMARK_E2F_TARGETS"               
# # [13] "HALLMARK_MYC_TARGETS_V1"             "HALLMARK_MYC_TARGETS_V2"            
# # [15] "HALLMARK_P53_PATHWAY"                "HALLMARK_KRAS_SIGNALING_UP"         
# # [17] "HALLMARK_KRAS_SIGNALING_DN"          "winlose_diff"                       
# # [19] "proapop_diff"                        "winapop_diff"                       
# # [21] "winloseapop_diff"    


# rc <- 2
# th <- 20
# complete <- FALSE

# res <- lapply(type, function(sel_type) {
#     print(sel_type)
#     selfiles <- list.files("res/res_score/", full.names = TRUE, pattern = str_glue("{sel_type}_E.*_{complete}_{rc}_{th}.rds"))
#     pbmcapply::pbmclapply(selfiles, function(s) {
#         readRDS(s)
#     }, mc.cores = 27) %>% rbindlist()
# })
# res %>% ll
# # 21

# for (i in 2:length(res)) {
#     if (i == 2) {
#         print(i)
#         new_res <- full_join(
#             select(res[[i - 1]], cellname, !contains("new")),
#             select(res[[i]], cellname, !contains("new"))
#         )
#     } else {
#         print(i)
#         new_res <- full_join(new_res, select(res[[i]], cellname, !contains("new")))
#     }
# }

# save_new_res <- new_res %>%
#     separate(cellname, c("corrd", "timepoint", "section"), sep = "-", remove = FALSE, convert = TRUE) %>%
#     separate(corrd, c("y", "x"), sep = "_", convert = TRUE) %>%
#     mutate(time = as.numeric(sub("E", "", timepoint))) %>%
#     rename(celltype2 = annotation) %>%
#     rename(annotation = celltype) %>%
#     rename(celltype = celltype2) %>% 
#     select(cellname, y, x, timepoint, time, section, celltype, annotation, n_grid, everything())
# saveRDS(save_new_res, "res/summary_entropy.rds")

# f_save_new_res <- save_new_res %>%
#     select(
#         -proliferation_moscot_entropy, -proliferation_moscot_genescore,
#         -essential_all_entropy, -essential_all_genescore,
#         -essential_both_entropy, -essential_both_genescore,
#         -haploinsufficiency_entropy, -haploinsufficiency_genescore,
#         -winloseapop_diff_entropy, -winloseapop_diff_genescore,
#         -winapop_diff_entropy, -winapop_diff_genescore,
#         -winloseapop_diff_entropy, -winloseapop_diff_genescore,
#         -HALLMARK_KRAS_SIGNALING_DN_entropy, -HALLMARK_KRAS_SIGNALING_DN_genescore,
#         -HALLMARK_KRAS_SIGNALING_UP_entropy, -HALLMARK_KRAS_SIGNALING_UP_genescore,
#         -HALLMARK_IL6_JAK_STAT3_SIGNALING_entropy, -HALLMARK_IL6_JAK_STAT3_SIGNALING_genescore,
#         -HALLMARK_MTORC1_SIGNALING_entropy, -HALLMARK_MTORC1_SIGNALING_genescore,
#         -HALLMARK_E2F_TARGETS_entropy, -HALLMARK_E2F_TARGETS_genescore, 
#         -HALLMARK_WNT_BETA_CATENIN_SIGNALING_entropy, -HALLMARK_WNT_BETA_CATENIN_SIGNALING_genescore,
#     )
# saveRDS(f_save_new_res, "res/f_summary_entropy.rds")

#######################################################################
### check scores (genescore, entropy, new values) #####################
#######################################################################

res <- readRDS("res/f_summary_entropy.rds")

# timescore data (gene score)
res_genescore <- res %>%
    select(time, ends_with("genescore")) %>%
    gather(key = "type", value = "genescore", -time) %>%
    mutate(type = factor(type)) %>%
    as_tibble()

res_genescore$type <- factor(res_genescore$type,
    levels = c(
        "win_genescore", "lose_genescore",
        "proliferation_plos_genescore", "HALLMARK_MYC_TARGETS_V1_genescore", "HALLMARK_MYC_TARGETS_V2_genescore",
        "apoptosis_msigdbhallmark_genescore", "HALLMARK_P53_PATHWAY_genescore",
        "winlose_diff_genescore",
        "proapop_diff_genescore"
    ),
    labels = c(
        "win", "lose",
        "proliferation", "MYC targets V1 (HALLMARK)", "MYC targets V2 (HALLMARK)",
        "apoptosis", "P53 pathway (HALLMARK)", 
        "competition rate",
        "growth rate"
    )
)


# ggplot(res_genescore, aes(time, genescore, color = type)) +
#     geom_smooth(method = "loess", se = FALSE) +
#     scale_x_continuous(breaks = seq(9.5, 16.5)) +
#     geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
#     scale_color_manual(values = col) +
#     theme_bw() +
#     theme(
#         axis.text = element_text(color = "black"),
#         strip.background = element_blank(),
#         strip.text = element_text(hjust = 0)
#     )
# ggsave("fig/genescore_time.png", units = "in", height = 5, width = 10)

medina_res_genescore <- res_genescore %>%
    group_by(time, type) %>%
    summarise(median_genescore = median(genescore, na.rm = TRUE)) %>%
    mutate(linetype = if_else(grepl("rate", type), "ori", "diff"))
col <- WGCNA::labels2colors(levels(medina_res_genescore$type))


ggplot(medina_res_genescore, aes(time, median_genescore, color = type, linetype = linetype)) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    # geom_smooth(method = "loess", se = FALSE) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = seq(9.5, 16.5)) +
    scale_color_manual(values = col) +
    theme_bw() +
    labs(y = "genescore (median)", x = "time point") +
    theme(
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
    ) +
    guides(linetype = "none")
ggsave("fig/genescore_time_median.png", units = "in", height = 4, width = 6)


# timescore data (entropy)
res_entropy <- res %>%
    select(time, ends_with("entropy")) %>%
    gather(key = "type", value = "entropy", -time) %>%
    mutate(type = factor(type)) %>%
    as_tibble()

# ggplot(res_entropy, aes(time, entropy, color = type)) +
#     geom_smooth(method = "loess", se = FALSE) +
#     scale_x_continuous(breaks = seq(9.5, 16.5)) +
#     geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
#     scale_color_manual(values = col) +
#     theme_bw() +
#     theme(
#         axis.text = element_text(color = "black"),
#         strip.background = element_blank(),
#         strip.text = element_text(hjust = 0)
#     )
# ggsave("fig/entropy_time.png", units = "in", height = 5, width = 10)


res_entropy$type <- factor(res_entropy$type,
    levels = c(
        "win_entropy", "lose_entropy",
        "proliferation_plos_entropy", "HALLMARK_MYC_TARGETS_V1_entropy", "HALLMARK_MYC_TARGETS_V2_entropy",
        "apoptosis_msigdbhallmark_entropy", "HALLMARK_P53_PATHWAY_entropy",
        "winlose_diff_entropy",
        "proapop_diff_entropy"
    ),
    labels = c(
        "win", "lose",
        "proliferation", "MYC targets V1 (HALLMARK)", "MYC targets V2 (HALLMARK)",
        "apoptosis", "P53 pathway (HALLMARK)",
        "competition rate",
        "growth rate"
    )
)

medina_res_entropy <- res_entropy %>%
    group_by(time, type) %>%
    summarise(median_entropy = median(entropy, na.rm = TRUE)) %>%
    mutate(linetype = if_else(grepl("rate", type), "ori", "diff"))

medina_res_entropy[apply(is.na(medina_res_entropy), 1, sum) != 0, ] %>% as.data.frame
# col <- WGCNA::labels2colors(levels(medina_res_entropy$type))

ggplot(medina_res_entropy, aes(time, median_entropy, color = type, linetype = linetype)) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    # geom_smooth(method = "loess", se = FALSE) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = seq(9.5, 16.5)) +
    scale_color_manual(values = col) +
    theme_bw() +
    labs(y = "entropy (median)", x = "time point") +
    theme(
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)
    ) +
    guides(linetype = "none")
ggsave("fig/entropy_time_median.png", units = "in", height = 4, width = 6)


#################################################################################
### Cor of Marker scores ########################################################
#################################################################################

rm(list = ls())
library(tidyverse)
library(data.table)

library(jjutil)
library(tidytext)


res <- readRDS("res/f_summary_entropy.rds")
colnames(res)[10:27]
colnames(res)[10:27] <- c(
    "win (entropy)", "win (genescore)",
    "lose (entropy)", "lose (genescore)",
    "proliferation  (entropy)", "proliferation (genescore)",
    "apoptosis (entropy)", "apoptosis (genescore)",
    "MYC targets V1 (HALLMARK) (entropy)", "MYC targets V1 (HALLMARK) (genescore)",
    "MYC targets V2 (HALLMARK) (entropy)", "MYC targets V2 (HALLMARK) (genescore)",
    "P53 pathway (HALLMARK) (entropy)", "P53 pathway (HALLMARK) (genescore)",
    "competition rate (entropy)", "competition rate(genescore)",
    "growth rate (entropy)", "growth rate (genescore)"
)

# all tisses
res %>%
    select(-c(1:9)) %>%
    cor(use = "complete.obs") %>%
    # ggcorrplot::ggcorrplot(hc.order = T, type = "lower", lab = TRUE) +
    ggcorrplot::ggcorrplot(lab = TRUE) +
    # ggtitle(str_glue("{seltime}")) +
    theme(axis.text = element_text(color = "black"))
ggsave(str_glue("fig/cor.png"), width = 12, height = 12, units = "in", bg = "white", limitsize = FALSE)



# seltime <- "E9.5"
# seltissue <- "Brain"
lapply(unique(res$timepoint), function(seltime) {  
    sel_res <- res %>%
        filter(timepoint == seltime)
    
    # all tisses
    sel_res %>%
        select(-c(1:9)) %>%
        cor(use = "complete.obs") %>%
        # ggcorrplot::ggcorrplot(hc.order = T, type = "lower", lab = TRUE) +
        ggcorrplot::ggcorrplot(lab = TRUE) +
        ggtitle(str_glue("{seltime}")) +
        theme(axis.text = element_text(color = "black"))
    ggsave(str_glue("fig/cor/cor_{seltime}.png"), width = 12, height = 12, units = "in", bg = "white", limitsize = FALSE)

    # separate tissue
    pbmcapply::pbmclapply(unique(sel_res$celltype), function(seltissue) {
        #cat(str_glue("{seltime}_{seltissue}"), "\n")
        sel_res %>%
            filter(celltype == seltissue) %>% 
            select(-c(1:9)) %>%
            cor(use = "complete.obs") %>% 
            # ggcorrplot::ggcorrplot(hc.order = T, type = "lower", lab = TRUE) + 
            ggcorrplot::ggcorrplot(lab = TRUE) +
            ggtitle(str_glue("{seltime}_{seltissue}")) +
            theme(axis.text = element_text(color = "black"))
        ggsave(str_glue("fig/cor/cor_{seltime}_{seltissue}.png"), width = 12, height = 12, units = "in", bg = "white", limitsize = FALSE)
    }, mc.cores = 25)
})
