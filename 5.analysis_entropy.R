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
    make_option("--th",
        action = "store", default = 10, type = "integer",
        help = "th for cal entropy"
    ),
    make_option("--numbin",
        action = "store", default = 10, type = "integer",
        help = "numbin for cal entropy"
    ),
    make_option("--class",
        action = "store", default = NULL, type = "character",
        help = "homo? hetero?"
    )
)

opt <- parse_args(OptionParser(option_list = option_list))

rc <- opt$rc
th <- opt$th
numbin <- opt$numbin

class <- opt$class

# rc <- 2
# th <- 10
# numbin <- 10

# class <- "hetero"
# sel <- "win"


# read entropy results
rawres <- readRDS(str_glue("res/summary_entropy_{rc}_{th}_{numbin}.rds"))

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

    sel_names <- str_glue(rc, "_", th, "_", numbin, "_", sel)
    cat(sel_names, "\n")
        
    # homo? or hetero?
    if (class == "homo") {
        sel_res <- res %>%
            filter(!grepl(";", celltype))
    } else {
        sel_res <- res %>%
            filter(grepl(";", celltype))
    }

    summary_res <- append(summary_res, list(sel_res))

    # boxplot color
    if (sel %in% c("ratio", "winlose_diff", "winapop_diff", "winloseapop_diff", "win", "lose")) {
        set_col <- c("white", "firebrick3")
    } else {
        set_col <- c("white", "royalblue3")
    }

    ggplot(sel_res, aes(entropy, reorder_within(celltype, entropy, timepoint, median), fill = col)) +
        theme_bw() +
        geom_boxplot() +
        scale_fill_manual(values = set_col) +
        facet_wrap(~timepoint, scales = "free", nrow = 4) +
        theme(legend.position = "none")

    if (class == "homo") {
        ggsave(str_glue("fig/entropy/entropy_{sel}_{rc}_{numbin}_{class}_{th}.png"), width = 10, height = 20, units = "in", bg = "white")
    } else {
        ggsave(str_glue("fig/entropy/entropy_{sel}_{rc}_{numbin}_{class}_{th}.png"), width = 15, height = 40, units = "in", bg = "white")
    }
}

saveRDS(rbindlist(summary_res), str_glue("data/entropy/entropy_{rc}_{numbin}_{class}_{th}.rds"))
