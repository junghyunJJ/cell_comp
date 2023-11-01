rm(list = ls())
library(tidyverse)
library(data.table)

library(jjutil)
library(Seurat)

#######################################################################
### resutls summary ###################################################
#######################################################################

# allfiles <- list.files("data/renormalized", pattern = "csv")

# raw_summary_mosta <- pbmcapply::pbmclapply(allfiles, function(sel) {
#     cat(sel, "\n")
#     fread(str_glue("data/renormalized/{sel}"))
# }, mc.cores = 25) %>% rbindlist

# summary_mosta <- raw_summary_mosta %>%
#     mutate(timepoint = factor(timepoint, levels = c("E9.5", "E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16.5"))) %>%
#     mutate(breaks_ratio = factor(breaks_ratio)) %>% 
#     mutate(annotation = factor(annotation)) %>% 
#     mutate(breaks_ratio = cut(ratio, breaks = 7, include.lowest = TRUE)) %>% 
#     arrange(time) %>% as_tibble
   
# # The "cell_name" column should use x and y cordinate.
# # NOTE!! There is a some extra naming in E15.5 (e.g., "_HC21_E15.5_main"), so remove.
# # cell_name = "newy"_"newx"-"timepoint"-"section"
# summary_mosta <- summary_mosta %>%
#     mutate(cell_name = sub("_HC21_E15.5_main", "", cell_name)) %>%
#     separate(cell_name, c("newy", "newx"), sep = "_", convert = TRUE, remove = FALSE) %>%
#     mutate(cell_name = str_glue("{newy}_{newx}-{timepoint}-{section}"))    

# saveRDS(summary_mosta, "res/summary_mosta.rds")

#######################################################################
### check data ########################################################
#######################################################################

summary_mosta <- readRDS("res/summary_mosta.rds")

summary_mosta %>%
    select("timepoint", "section") %>%
    unique.data.frame() %>%
    group_by(timepoint) %>%
    count()
# timepoint `sum(n)`
# <fct>        <int>
# 1 E9.5             5
# 2 E10.5            4
# 3 E11.5            4
# 4 E12.5            6
# 5 E13.5            4
# 6 E14.5            7
# 7 E15.5            5
# 8 E16.5           18

summary_mosta %>%
    select("timepoint", "section") %>%
    unique.data.frame() %>%
    group_by(timepoint) %>%
    count() %>%
    ungroup() %>%
    summarise(sum(n))
# a total 53 section across 8 timepoints

summary_mosta$annotation %>% unique %>% ll
# 49 annotation


# #######################################################################
# ### check coord and annotation ########################################
# #######################################################################
sel <- "E9.5"
for (sel in levels(summary_mosta$timepoint)) {
    cat(sel, "\n")
    ggori <- summary_mosta %>% filter(timepoint == sel) %>%
        ggplot(aes(x, y, col = annotation)) +
        geom_point(size = 0.5) +
        theme_bw() +
        facet_wrap(~section, scales = "free", nrow = 1)

  
    ggnew <- summary_mosta %>% filter(timepoint == sel) %>% 
        ggplot(aes(newx, newy, col = annotation)) +
        geom_point(size = 0.5) +
        theme_bw() +
        facet_wrap(~section, scales = "free", nrow = 1)
  
    cowplot::plot_grid(ggori, ggnew, nrow = 2)

    if (sel == "E16.5") {
        ggsave(filename = paste0("fig/check/", sel, ".png"), width = 60, height = 10, units = "in", limitsize = FALSE)    
    } else if (sel == "E9.5") {
        ggsave(filename = paste0("fig/check/", sel, ".png"), width = 20, height = 10, units = "in", limitsize = FALSE)    
    }else {
        ggsave(filename = paste0("fig/check/", sel, ".png"), width = 30, height = 10, units = "in", limitsize = FALSE)
    }
}
