rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)

library(jjutil)
library(Seurat)


allfiles_mh <- list.files("data/renormalized", pattern = "_obs", full.names = TRUE)

# NOTE!!!! some files have little bit different "cell_name"
res_mh <- pbmcapply::pbmclapply(allfiles_mh, function(d) {
    tmp_d <- fread(d)
    colnames(tmp_d)[1] <- "cell_name"
    tmp_d %>%
        mutate(cell_name = sub("([0-9]+_[0-9]+).*", "\\1", cell_name)) %>%
        mutate(timepoint = sub(".*(E.*)_(E.*).MOSTA_renormalized_obs.csv", "\\1", d)) %>%
        mutate(section = sub(".*(E.*)_(E.*).MOSTA_renormalized_obs.csv", "\\2", d)) %>%
        mutate(cell_name = paste0(cell_name, "-", timepoint, "-", section))
}, mc.cores = 20) %>% rbindlist()

# we only focused on the hallmarkdata 
f_res_mh <- res_mh %>% select("cell_name", starts_with("HALLMARK_"))
f_res_mh %>% h

# the only differnet between "summary_mosta.rds" and "final_summary_mosta.rds" is that "newy" and "newx" column
# summary_mosta <- readRDS("res/summary_mosta.rds")
summary_mosta <- readRDS("res/final_summary_mosta.rds")
summary_mosta <- summary_mosta %>%
    mutate(cell_name = sub("([0-9]+_[0-9]+)\\.(E.*)", "\\1-\\2", cell_name)) %>%
    mutate(cell_name = paste0(cell_name, "-", section))
summary_mosta %>% h

# NOTE!!!!!!
# cell_name = "newy"_"newx"-"timepoint"-"section"
dat <- inner_join(summary_mosta, f_res_mh)

# PERFCET!!
nrow(f_res_mh)
nrow(summary_mosta)
nrow(dat)

saveRDS(dat, "res/new_summary_mosta.rds")