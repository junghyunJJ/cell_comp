rm(list = ls())
library(tidyverse)
library(data.table)

library(jjutil)
library(optparse)

source("functions.R")

option_list <- list(
    make_option("--rc",
        action = "store", default = 2, type = "integer",
        help = "window width for calculate entropy"
    ),
    make_option("--numbin",
        action = "store", default = 10, type = "integer",
        help = "numbin for cal entropy"
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
#sel <- 1; sel_timepoint <- "E9.5"; sel_type <- "proliferation_plos";z <- 1;sel_section <- "E1S1";rc<-2;numbin<-10

sel_mosta <- mosta %>% 
    filter(timepoint == sel_timepoint)

section <- unique(sel_mosta$section)
#z <- 1;sel_section <- section[z]

res <- lapply(section, function(sel_section) {
    cat("type:", sel_type, "/ timepoint:", sel_timepoint, "/ section:", sel_section, "\n")

    # reorganization of data fotmat 
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

    # matrix data based on the x and y coordinate
    matdat <- make_mat(dat)
    
    # cal entropy based on the specific window (default = 2)
    s_matdat <- split_mat(matdat, r = opt$rc, c = opt$rc)
    celltype_comb <- make_celltype_comb(unique(dat$annotation), includehomo = TRUE)        
    res_entropy <- cal_entropy(s_matdat, celltype_comb, numbin = opt$numbin, unit = "log2")
    
    # summary of data
    res_entropy %>%
        mutate(timepoint = get("sel_timepoint")) %>%
        mutate(section = get("sel_section")) %>%
        mutate(type = get("sel_type"))
    
}) %>% rbindlist()

saveRDS(res, str_glue("res/res_entropy/{sel_type}_{sel_timepoint}_{opt$rc}_{opt$numbin}.rds"))
