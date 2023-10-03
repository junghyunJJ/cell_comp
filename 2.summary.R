rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)

library(jjutil)
library(Seurat)


# cal_ratio function from Patrick
cal_ratio <- function(norm_counts, win, lose) {
    weight <- rep(1, length(c(win, lose)))
    names(weight) <- c(win, lose)
    
    # weight for more important genes
    weight[names(weight) %in% c("Myc", "Mycn", "Ras", "Egfr")] <- 2
    weight[names(weight) == "Cacfd"] <- 3

    win_counts <- norm_counts[rownames(norm_counts) %in% win, ]
    win_counts %>% dim
    lose_counts <- norm_counts[rownames(norm_counts) %in% lose, ]
    lose_counts %>% dim

    for (i in seq_along(weight)) {
        if (names(weight)[i] %in% rownames(win_counts)) {
            win_counts[rownames(win_counts) == names(weight)[i], ] <- win_counts[rownames(win_counts) == names(weight)[i], ] * weight[i]
        } else if (names(weight)[i] %in% rownames(lose_counts)) {
            lose_counts[rownames(lose_counts) == names(weight)[i], ] <- lose_counts[rownames(lose_counts) == names(weight)[i], ] * weight[i]
        } else {
            next
        }
    }

    win_counts <- colMeans(win_counts)
    lose_counts <- colMeans(lose_counts)

    ratio <- (win_counts + 1) / (lose_counts + 1) - 1 # NOTE!!!! any reference?? or related papers? -> NO!
    breaks <- cut(ratio, breaks = 7, include.lowest = TRUE)
    return(list(ratio = ratio, breaks = breaks))
}

#############################################################################
#############################################################################
#############################################################################

allfiles <- list.files("data/renormalized", pattern = "MOSTA")
# allfiles <- allfiles[1:2]

win <- readLines("data/win_mouse.txt")
lose <- readLines("data/lose_mouse.txt")

summary_mosta <- lapply(allfiles, function(sel) {
    
    #sel <- "E9.5_E1S1.MOSTA_renormalized.h5ad"
    cat(sel, "\n")
        
    # There is a erro when we use LoadH5Seurat function in SeuratDisk (https://github.com/mojaveazure/seurat-disk/issues/109)
    # So, we directly read the scanpy output (i.e., h5ad) using anndata::read_h5ad.
    dat <- anndata::read_h5ad(str_glue("data/renormalized/{sel}"))


    # 1. exp data 
    norm_counts <- dat$X %>% as.matrix %>% t
    res_cal_ratio <- cal_ratio(norm_counts, win, lose)


    # 2. coord data
    coord <- dat$obsm$spatial %>% as.data.frame %>% 
        set_colnames(c("x", "y")) %>% 
        #mutate(x = as.integer(x)) %>% 
        #mutate(y = as.integer(abs(y))) %>% 
        mutate(cell_name = rownames(dat$obs)) %>% 
        select(cell_name, everything()) %>% as_tibble() %>% 

        mutate(ratio = res_cal_ratio$ratio) %>% 
        mutate(breaks_ratio = res_cal_ratio$breaks)


    # 3. meta data
    # set names
    tmp <- unlist(strsplit(sel, "_"))
    timepoint <- tmp[1]
    time <- sub("E", "", timepoint) %>% as.numeric
    section <- sub(".MOSTA", "", tmp[2])

    meta <- dat$obs %>% 
        select(annotation, starts_with(c("win", "lose", "proliferation", "apoptosis", "essential", "haploinsufficiency"))) %>% 
        rownames_to_column("cell_name") %>% 
        mutate(timepoint = timepoint) %>% 
        mutate(time = time) %>%
        mutate(section = section) %>%
        select(cell_name, timepoint, time, section, everything()) %>% as_tibble()

    res <- inner_join(coord, meta, by = "cell_name")
    res
}) %>% rbindlist()

summary_mosta <- summary_mosta %>% 
    mutate(timepoint = factor(timepoint, levels = c("E9.5", "E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16.5"))) %>% 
    mutate(breaks_ratio = factor(breaks_ratio)) %>% 
    mutate(annotation = factor(annotation)) %>% 
    mutate(breaks_ratio = cut(ratio, breaks = 7, include.lowest = TRUE)) %>% 
    arrange(time) %>% as_tibble

saveRDS(summary_mosta, "res/summary_mosta.rds")