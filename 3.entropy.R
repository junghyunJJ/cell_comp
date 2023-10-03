rm(list = ls())
library(tidyverse)
library(data.table)

library(jjutil)
library(optparse)

source("functions.R")

option_list <- list(
    make_option("--analysis", action = "store", default = FALSE, type = "character", 
        help = "please select type of anlaysis 1. window (please also set r and c for window size) or 2. tissue (whole tissue) or 3. edege"
    ),
    make_option("--rc", action = "store", default = NA, type = "integer", 
        help = "window width for calculate entropy"
    )
)

opt <- parse_args(OptionParser(option_list = option_list))

rawmosta <- readRDS("res/new_summary_mosta.rds")

scaling <- 1
mosta <- rawmosta %>%
    mutate(winlose_diff = (win - lose) / scaling) %>%
    mutate(proapop_diff = (proliferation_plos - apoptosis_msigdbhallmark) / scaling) %>%
    mutate(winapop_diff = (win - apoptosis_msigdbhallmark) / scaling) %>% 
    mutate(winloseapop_diff = (win - apoptosis_msigdbhallmark - lose) / scaling)

# set for array in Slurm
all_timepoint <- mosta$timepoint %>% unique %>% as.character()
all_type <- colnames(mosta)[c(6, 12:ncol(mosta))]
(all <- expand.grid(all_timepoint, all_type))
# 512

#######################################################################################
#######################################################################################
#######################################################################################
sel <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
(sel_timepoint <- all[sel, 1] %>% as.character())
(sel_type <- all[sel, 2] %>% as.character())
# sel <- 1; sel_timepoint <- "E9.5"; sel_type <- "proliferation_plos";z <- 1;sel_section <- "E1S1";c<-r<-2

sel_mosta <- mosta %>% 
    filter(timepoint == sel_timepoint)

section <- unique(sel_mosta$section)
# z <- 1;sel_section <- section[z]

res <- lapply(section, function(sel_section) {
    cat("type:", sel_type, "/ timepoint:", sel_timepoint, "/ section:", sel_section, "\n")

    # data fotmatreorganization
    dat <- sel_mosta %>%
        filter(section == get("sel_section")) %>%
        select(newx, newy, annotation, get("sel_type"), timepoint, section) %>%
        mutate(x = newx + 1) %>%
        mutate(y = newy + 1) %>%
        rename(value = get("sel_type")) %>%
        select(x, y, annotation, value, everything()) %>%
        arrange(x, y)
    
    #  transformation of cell score data using "trans_dat" function (the value range is 0-1)
    dat <- dat %>%
        rename(ori_value = value) %>%
        mutate(value = trans_dat(ori_value))

    matdat <- make_mat(dat)
    
    # type of analysis
    if (opt$analysis == "window") {
        r <- c <- opt$rc
        # cal entropy based on the specific window (or bin)
        s_matdat <- split_mat(matdat, r = r, c = c)
        rawres <- cal_entropy(matdat$mat_dat, s_matdat, r = r, c = c, numbin = 5)

    } else if (opt$analysis == "tissue") {

        celltype <- unique(dat$annotation)
        rawres <- lapply(celltype, function(sel_celltype) {
            raw_sel_dat <- matdat$mat_dat[matdat$anno_dat == sel_celltype]
            sel_dat <- raw_sel_dat[!is.na(raw_sel_dat)]
            entropy::discretize(sel_dat, numBins = 5, r = c(0, 1)) %>%
                entropy::entropy.empirical(unit = "log2")
        }) %>% unlist
        rawres <- data.frame(celltype = celltype, entropy = rawres)
    }

    # summary of data
    rawres %>%
        # filter(!grepl(";", annotation)) %>%
        mutate(timepoint = get("sel_timepoint")) %>%
        mutate(section = get("sel_section")) %>%
        mutate(type = get("sel_type"))
    
}) %>% rbindlist()

if (opt$analysis == "window") {
    saveRDS(res, str_glue("res/res_entropy/{sel_type}_{sel_timepoint}_{opt$analysis}_{opt$rc}.rds"))
} else if (opt$analysis == "tissue") {
    saveRDS(res, str_glue("res/res_entropy/{sel_type}_{sel_timepoint}_{opt$analysis}.rds"))
}