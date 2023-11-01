rm(list = ls())
library(tidyverse)
library(data.table)
library(tidytext)

library(jjutil)
library(openxlsx)

# library(Vennerable)
library(ggvenn)
library(fgsea)

res <- readRDS("res/final_summary_entropy.rds")
res %>% colnames()
load("res/csea_genescore.Rdata")

list_res_cellgroup[[4]]
sel_cell <- list_res_cell[[4]]

sel_dat <- res %>%
    filter(cellname %in% sel_cell) %>%
    filter(grepl("Heart", annotation))
unique(sel_dat$section) %>% unique()

sel_dat$cellname
sel_timepoint <- "E12.5"
#####################################################
### plot ############################################
#####################################################
# summary_mosta <- readRDS("res/summary_mosta.rds")

# plot_sel_data <- filter(summary_mosta, cell_name %in% sel_dat$cellname)

# summary_mosta %>%
#     filter(section == unique(sel_dat$section), timepoint == "E12.5") %>%
#     ggplot(aes(x, -y, col = annotation)) +
#     geom_point(size = 0.5, alpha = 0.2) +
#     geom_point(data = plot_sel_data, aes(x, -y), col = "red") +
#     theme_void()
# ggsave("fig/E12.5_sel.png", height = 8, width = 8, bg = "white")

# summary_mosta %>%
#     filter(section == unique(sel_dat$section), timepoint == "E12.5") %>%
#     ggplot(aes(x, -y, col = annotation)) +
#     geom_point(size = 0.5) +
#     # geom_point(size = 0.5, alpha = 0.2) +
#     # geom_point(data = plot_sel_data, aes(x, -y), col = "red") +
#     theme_void()
# ggsave("fig/ori_E12.5_sel.png", height = 8, width = 8, bg = "white")


#####################################################
#####################################################
#####################################################
sel_dat %>% dim

i <- 1

res_cs <- lapply(seq_len(nrow(sel_dat)), function(i) {
    # select a center cell (c_cells)
    c_dat <- sel_dat[i, ]
    x <- c_dat$x
    y <- c_dat$y
    
    # select surrounding cells (s_cells)
    s_coord <- list(
        c(y + 1, x), # 12'
        c(y + 1, x + 1), # 1.5'
        c(y, x + 1), # 3'
        c(y - 1, x + 1), # 4.5'
        c(y - 1, x), # 6'
        c(y - 1, x - 1), # 7.5'
        c(y, x - 1), # 9'
        c(y + 1, x - 1) # 10.5'
    )
    s_coord <- lapply(s_coord, paste0, collapse = "_") %>% unlist()
    s_dat <- res %>%
        filter(section == unique(sel_dat$section), timepoint == "E12.5") %>%
        unite("coord", y:x, sep = "_", remove = FALSE) %>%
        filter(coord %in% s_coord) %>%
        select(-coord, -x, -y, -annotation, -n_grid)
    # browser()

    # check celltype: we only focued on the homo
    check_celltype <- c(c_dat$celltype, s_dat$celltype) %>% unique()
    
    # we skip res if there are NA in >50%
    # if (length(check_celltype) == 1 & sum(!is.na(s_dat$winlose_diff_entropy)) > 3 & sum(!is.na(s_dat$winlose_diff_genescore)) > 3) {
    # median value of surrounding cells
    median_s_dat <- s_dat %>%
        select(ends_with("entropy"), ends_with("genescore")) %>%
        apply(2, median, na.rm = TRUE)
    
    # median_surrounding_dat - center_dat
    value_c_dat <- c_dat %>% select(ends_with("entropy"), ends_with("genescore"))
    delta_dat <- median_s_dat - value_c_dat
    colnames(delta_dat) <- str_glue("delta_{colnames(delta_dat)}")

    summary_dat <- cbind(c_dat[, c(1:4, 6, 7)],
        delta_dat,
        n_s_cells_genescore = sum(!is.na(s_dat$winlose_diff_genescore)),
        n_s_cells_entropy = sum(!is.na(s_dat$winlose_diff_entropy))
    )
    # } else {
    #     summary_dat <- NA
    # }
    summary_dat
}) %>% rbindlist()

t_res_cs <- res_cs %>%
    select(cellname, celltype, ends_with("entropy"), ends_with("genescore")) %>% 
    gather(key = "type", value = "value", -cellname, -celltype, -n_s_cells_genescore, -n_s_cells_entropy) %>%
    mutate(score = if_else(grepl("entropy", type), "entropy", "genescore")) %>%
    mutate(type = sub("(_entropy|_genescore)", "", type))

t_res_cs$type <- factor(t_res_cs$type,
    levels = c(
        "delta_win", "delta_lose",
        "delta_proliferation_plos", "delta_HALLMARK_MYC_TARGETS_V1", "delta_HALLMARK_MYC_TARGETS_V2",
        "delta_apoptosis_msigdbhallmark", "delta_HALLMARK_P53_PATHWAY",
        "delta_winlose_diff",
        "delta_proapop_diff"
    ) %>% rev(),
    labels = c(
        "win", "lose",
        "proliferation", "MYC targets V1 (HALLMARK)", "MYC targets V2 (HALLMARK)",
        "apoptosis", "P53 pathway (HALLMARK)",
        "competition rate",
        "growth rate"
    ) %>% rev()
)

    
ggplot(t_res_cs, aes(value, type)) +
    geom_boxplot() +
    theme_bw() +
    scale_fill_manual(values = c("grey60", "firebrick3")) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    facet_grid(celltype ~ score, scales = "free_x") +
    labs(y = "", x = "expression level of\nmedian value of surrounding cells (8 cells) - center (1 cell)") +
    theme(
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)
    )
ggsave(str_glue("fig/new_apoptosis/{sel_timepoint}_apoptosis.png"), width = 8, height = 4, unit = "in", limitsize = FALSE)
å