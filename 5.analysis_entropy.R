rm(list = ls())
library(tidyverse)
library(data.table)
library(tidytext)

library(jjutil)
library(openxlsx)

library(optparse)

option_list <- list(
    make_option("--rc",
        action = "store", default = 2, type = "integer",
        help = "window width for calculate entropy"
    ),
    make_option("--numbin",
        action = "store", default = 10, type = "integer",
        help = "numbin for cal entropy"
    ),
    make_option("--class",
        action = "store", default = NULL, type = "character",
        help = "homo? hetero?"
    ),
    make_option("--th",
        action = "store", default = NULL, type = "integer",
        help = "th for vislization"
    )
)

opt <- parse_args(OptionParser(option_list = option_list))

rc <- opt$rc
numbin <- opt$numbin
class <- opt$class
th <- opt$th

# rc <- 2
# numbin <- 10
# class <- "homo"

# sel <- "win"
# th <- 5

# read entropy results
rawres <- readRDS(str_glue("res/summary_entropy_2_{numbin}.rds"))

res <- rawres %>%
    filter(type %in% c(
        "ratio", "winlose_diff", "proapop_diff", "winapop_diff", "winloseapop_diff",
        "win", "lose", "proliferation_plos", "apoptosis_msigdbhallmark", "essential_both", "haploinsufficiency",
        "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"
    )) %>%
    mutate(timepoint = factor(timepoint, levels = c("E9.5", "E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16.5"))) %>%
    # mutate(col = if_else(celltype %in% c("Brain", "Spinal cord"), "1", "0")) %>%
    mutate(col = if_else(grepl("Brain|Spinal cord", celltype), "1", "0")) %>%
    mutate(type = factor(type, levels = c(
        "ratio", "winlose_diff", "proapop_diff", "winapop_diff", "winloseapop_diff",
        "win", "lose", "proliferation_plos", "apoptosis_msigdbhallmark", "essential_both", "haploinsufficiency",
        "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"
    ), 
    labels = c(
        "ratio", "winlose_diff", "proapop_diff", "winapop_diff", "winloseapop_diff",
        "win", "lose", "proliferation", "apoptosis", "essential", "haploinsufficiency",
        "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"
    ))) %>%
    mutate(timepoint = factor(timepoint, levels = c("E9.5", "E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16.5"))) %>% 
    arrange(timepoint)

##############################################################################
### plot #####################################################################
##############################################################################

summary_res <- list()

for (sel in levels(res$type)) {

    sel_names <- str_glue(rc, "_", numbin, "_", sel)
    cat(sel_names, "\n")

    sel_res <- res %>%
        filter(type == sel, !is.na(entropy)) %>%
        group_by(timepoint, type, celltype) %>%
        summarise(n = n(), median_entropy = median(entropy), .groups = "keep") %>%
        arrange(desc(n)) %>% 
        ungroup() %>%
        filter(n > th) # NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    # homo? or hetero?
    if (class == "homo") {
        sel_res <- sel_res %>%
            # filter(grepl(";", celltype)) %>%
            # filter(str_count(celltype, ";") == 1) %>%
            filter(!grepl(";", celltype)) %>%
            arrange(desc(median_entropy))
    } else {
        sel_res <- sel_res %>%
            filter(grepl(";", celltype)) %>%
            filter(str_count(celltype, ";") == 1) %>%
            arrange(desc(median_entropy))
    }

    # filtering based on the annotaiton of cell type
    f_sel_res <- res %>%
        filter(type == sel) %>%
        filter(celltype %in% sel_res$celltype) %>%
        mutate(th = th) %>%
        mutate(class = class)
    summary_res <- append(summary_res, list(f_sel_res))

    # boxplot color
    if (sel %in% c("ratio", "winlose_diff", "proapop_diff", "winapop_diff", "winloseapop_diff", "win", "lose")) {
        set_col <- c("white", "firebrick3")
    } else {
        set_col <- c("white", "royalblue3")
    }

    gg_sel_res <- ggplot(f_sel_res, aes(entropy, reorder_within(celltype, entropy, timepoint, median), fill = col)) +
        theme_bw() +
        geom_boxplot() +
        scale_fill_manual(values = set_col) +
        facet_wrap(~timepoint, scales = "free", nrow = 4) +
        theme(legend.position = "none")

    if (class == "homo") {
        ggsave(str_glue("fig/entropy/entropy_window_{rc}_{sel}_{class}_{th}.png"), width = 15, height = 20, units = "in", bg = "white")
    } else {
        ggsave(str_glue("fig/entropy/entropy_window_{rc}_{sel}_{class}_{th}.png"), width = 15, height = 40, units = "in", bg = "white")
    }
}

saveRDS(rbindlist(summary_res), str_glue("data/entropy/entropy_window_{rc}_{class}_{th}.rds"))
