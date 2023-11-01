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
    make_option("--th",
        action = "store", default = 20, type = "integer",
        help = "th for cal entropy"
    ),
    make_option("--complete",
        action = "store_true", default = FALSE,
        help = "only using complet data"
    )
)
opt <- parse_args(OptionParser(option_list = option_list))

mosta <- readRDS("res/summary_mosta.rds")

#######################################################################################
### set param for Slurm array #########################################################
#######################################################################################

# set param file for Slurm array
tmp_type <- grep("_entropy$", colnames(mosta), value = TRUE)
type <- sub("_entropy", "", tmp_type)

all_timepoint <- mosta$timepoint %>% unique %>% as.character()
tmp_param <- mosta %>% select(timepoint, section) %>% unique.data.frame()

param <- data.frame(
    timepoint = rep(tmp_param$timepoint, each = length(type)),
    section = rep(tmp_param$section, each = length(type)),
    type = rep(type, 53)
)
param %>% dim
# 1113

#######################################################################################
#######################################################################################
#######################################################################################

sel <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# sel <- 1;rc <- 2;sel_name<-"entropy";opt$complete<-FALSE
sel_timepoint <- param[sel, 1] %>% as.character()
sel_section <- param[sel, 2] %>% as.character()
sel_type <- param[sel, 3] %>% as.character()

cat("type:", sel_type, "/ section:", sel_section, "/ timepoint:", sel_timepoint, "/ complete:", opt$complete, "\n")

sel_mosta <- mosta %>%
    filter(timepoint == sel_timepoint, section == sel_section) %>% as_tibble
    
    
raw_res <- lapply(c("entropy", "genescore", "new", "new2"), function(sel_name) {
    sel_ana <- paste0(sel_type, "_", sel_name)

    # reorganization of data fotmat
    dat <- sel_mosta %>%
        dplyr::select(cell_name, newx, newy, annotation, any_of(sel_ana), timepoint, section) %>%
        rename(cellname = cell_name) %>% 
        mutate(x = newx + 1) %>%
        mutate(y = newy + 1) %>%
        rename(value = any_of(sel_ana)) %>% # please use "any_of" function
        select(x, y, annotation, value, everything()) %>%
        arrange(x, y)
    
    # matrix data based on the x and y coordinate
    matdat <- make_mat(dat)
  
    # cal entropy based on the specific window (default = 2)
    g_matdat <- grid_mat(matdat, r = opt$rc, c = opt$rc, complete = opt$complete)
    celltype_comb <- make_celltype_comb(g_matdat, th = opt$th, includehomo = TRUE, exclude = NULL)
    
    # summary of data
    res_celltype_comb <- summary_score(g_matdat, celltype_comb, name = sel_ana) %>%
        mutate(cellname = str_glue("{cellname}-{sel_timepoint}-{sel_section}")) %>%
        mutate(timepoint = sel_timepoint) %>%
        mutate(section = sel_section)        

    right_join(select(dat, cellname, annotation), res_celltype_comb, by = "cellname")
    
})
# lapply(raw_res, dim)
res <- cbind(raw_res[[1]], raw_res[[2]][, 4], raw_res[[3]][, 4], raw_res[[4]][, 4])
res <- res[apply(is.na(res), 1, sum) == 0, ]

saveRDS(res, str_glue("res/res_score/{sel_type}_{sel_timepoint}_{sel_section}_{opt$complete}_{opt$rc}_{opt$th}.rds"))
