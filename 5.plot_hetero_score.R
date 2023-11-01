rm(list = ls())
library(tidyverse)
library(data.table)

library(jjutil)
library(tidytext)

#######################################################################
### 1-0. prep genescore ###############################################
#######################################################################

res <- readRDS("res/f_summary_entropy.rds") 
res %>% colnames

# timescore data (gene score)
res_genescore <- res %>%
    select(time, annotation, ends_with("genescore")) %>%
    gather(key = "type", value = "genescore", -time, -annotation) %>%
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

res_genescore
res_genescore$type %>% unique

##############################################################################
### 1-1. genescore homo ######################################################
##############################################################################

# we only focus on the more than half of time points results (timepoints > 4)
sel_homo_tissue <- res_genescore %>%
    filter(!grepl(";", annotation)) %>%
    select(type, annotation, time) %>%
    unique.data.frame() %>%
    select(annotation, time) %>%
    table() %>%
    as.data.frame() %>%
    filter(Freq == 9) %>% # 9 types
    group_by(annotation) %>%
    count() %>%
    filter(n > 3) %>%
    select(annotation) %>%
    pull() %>%
    as.character()
sel_homo_tissue

homo_res_genescore <- res_genescore %>%
    filter(annotation %in% sel_homo_tissue)

median_homo_res_genescore <- homo_res_genescore %>%
    group_by(time, annotation, type) %>%
    summarise(median_genescore = median(genescore, na.rm = TRUE)) %>%
    mutate(linetype = if_else(grepl("rate", type), "ori", "diff"))
col <- WGCNA::labels2colors(levels(median_homo_res_genescore$type))

ggplot(median_homo_res_genescore, aes(time, median_genescore, color = type, linetype = linetype)) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = seq(9.5, 16.5)) +
    scale_color_manual(values = col) +
    facet_wrap(~ annotation, nrow = 6) +
    theme_bw() +
    labs(y = "genescore (median)", x = "time point") +
    theme(
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)
    ) +
    guides(linetype = "none")
ggsave("fig/homo_genescore_median.png", height = 15, width = 10, limitsize = FALSE)

##############################################################################
### 1-2. genescore hetero ####################################################
##############################################################################

# we only focus on the more than half of time points results (timepoints > 4)
sel_hetero_tissue <- res_genescore %>%
    filter(grepl(";", annotation)) %>%
    select(type, annotation, time) %>%
    unique.data.frame() %>%
    select(annotation, time) %>%
    table() %>%
    as.data.frame() %>%
    filter(Freq == 9) %>% # 9 types
    group_by(annotation) %>%
    count() %>%
    filter(n > 3) %>%
    select(annotation) %>%
    pull() %>%
    as.character()
sel_hetero_tissue

hetero_res_genescore <- res_genescore %>%
    filter(annotation %in% sel_hetero_tissue)


median_hetero_res_genescore <- hetero_res_genescore %>%
    group_by(time, annotation, type) %>%
    summarise(median_genescore = median(genescore, na.rm = TRUE)) %>%
    mutate(linetype = if_else(grepl("rate", type), "ori", "diff"))
# col <- WGCNA::labels2colors(levels(median_hetero_res_genescore$type))
# col[7] <- "red"
# col[8] <- "green"

ggplot(median_hetero_res_genescore, aes(time, median_genescore, color = type, linetype = linetype)) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = seq(9.5, 16.5)) +
    scale_color_manual(values = col) +
    facet_wrap(~ annotation, nrow = 6) +
    theme_bw() +
    labs(y = "genescore (median)", x = "time point") +
    theme(
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)
    ) +
    guides(linetype = "none")
ggsave("fig/hetero_genescore_median.png", height = 15, width = 20, limitsize = FALSE)

#######################################################################
### 2-0. prep entropy #################################################
#######################################################################

rm(list = ls())
res <- readRDS("res/f_summary_entropy.rds") 
res %>% colnames


# timescore data (entropy)
res_entropy <- res %>%
    select(time, annotation, ends_with("entropy")) %>%
    gather(key = "type", value = "entropy", -time, -annotation) %>%
    mutate(type = factor(type)) %>%
    as_tibble()

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
res_entropy
res_entropy$type %>% unique


##############################################################################
### 2-1. entropy homo ########################################################
##############################################################################

# we only focus on the more than half of time points results (timepoints > 4)
sel_homo_tissue <- res_entropy %>%
    filter(!grepl(";", annotation)) %>%
    select(type, annotation, time) %>%
    unique.data.frame() %>%
    select(annotation, time) %>%
    table() %>%
    as.data.frame() %>%
    filter(Freq == 9) %>% # 9 types
    group_by(annotation) %>%
    count() %>%
    filter(n > 3) %>%
    select(annotation) %>%
    pull() %>%
    as.character()
sel_homo_tissue

homo_res_entropy <- res_entropy %>%
    filter(annotation %in% sel_homo_tissue)

median_homo_res_entropy <- homo_res_entropy %>%
    group_by(time, annotation, type) %>%
    summarise(median_entropy = median(entropy, na.rm = TRUE)) %>%
    mutate(linetype = if_else(grepl("rate", type), "ori", "diff"))
col <- WGCNA::labels2colors(levels(median_homo_res_entropy$type))
# col[7] <- "red"
# col[8] <- "green"

ggplot(median_homo_res_entropy, aes(time, median_entropy, color = type, linetype = linetype)) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = seq(9.5, 16.5)) +
    scale_color_manual(values = col) +
    facet_wrap(~ annotation, nrow = 6) +
    theme_bw() +
    labs(y = "entropy (median)", x = "time point") +
    theme(
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)
    ) +
    guides(linetype = "none")
ggsave("fig/homo_entropy_median.png", height = 15, width = 10, limitsize = FALSE)

##############################################################################
### 2-2. entropy hetero ######################################################
##############################################################################

# we only focus on the more than half of time points results (timepoints > 4)
sel_hetero_tissue <- res_entropy %>%
    filter(grepl(";", annotation)) %>%
    select(type, annotation, time) %>%
    unique.data.frame() %>%
    select(annotation, time) %>%
    table() %>%
    as.data.frame() %>%
    filter(Freq == 9) %>% # 9 types
    group_by(annotation) %>%
    count() %>%
    filter(n > 3) %>%
    select(annotation) %>%
    pull() %>%
    as.character()
sel_hetero_tissue

hetero_res_entropy <- res_entropy %>%
    filter(annotation %in% sel_hetero_tissue)


median_hetero_res_entropy <- hetero_res_entropy %>%
    group_by(time, annotation, type) %>%
    summarise(median_entropy = median(entropy, na.rm = TRUE)) %>%
    mutate(linetype = if_else(grepl("rate", type), "ori", "diff"))
# col <- WGCNA::labels2colors(levels(median_hetero_res_entropy$type))
# col[7] <- "red"
# col[8] <- "green"

ggplot(median_hetero_res_entropy, aes(time, median_entropy, color = type, linetype = linetype)) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = seq(9.5, 16.5)) +
    scale_color_manual(values = col) +
    facet_wrap(~ annotation, nrow = 6) +
    theme_bw() +
    labs(y = "entropy (median)", x = "time point") +
    theme(
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)
    ) +
    guides(linetype = "none")
ggsave("fig/hetero_entropy_median.png", height = 15, width = 20, limitsize = FALSE)
