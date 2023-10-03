rm(list = ls())
library(tidyverse)
library(data.table)

library(jjutil)

# The "cell_name" column should use x and y cordinate.
# NOTE!! There is a some extra naming in E15.5 (e.g., "_HC21_E15.5_main"), so remove.

summary_mosta <- readRDS("res/summary_mosta.rds") %>% 
    mutate(cell_name = sub("_HC21_E15.5_main", "", cell_name)) %>% 
    separate(cell_name, c("newy", "newx"), sep = "_", convert = TURE, remove = FALSE) %>% 
    mutate(cell_name = paste0(cell_name, ".", timepoint))
saveRDS(summary_mosta, "res/final_summary_mosta.rds")


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



mosta %>% select("timepoint","section") %>% unique.data.frame() %>% group_by(timepoint) %>% count()
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

mosta %>% select("timepoint","section") %>% unique.data.frame() %>% group_by(timepoint) %>% count() %>% ungroup %>% summarise(sum(n)) 
# a total 53 section across 8 timepoints

mosta$annotation %>% unique %>% ll
# 49 annotation