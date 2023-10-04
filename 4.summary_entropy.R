rm(list = ls())
library(tidyverse)
library(data.table)

library(jjutil)
library(tidytext)

# scaling <- 1
# rawmosta <- readRDS("res/new_summary_mosta.rds")
# mosta <- rawmosta %>%
#     mutate(winlose_diff = (win - lose) / scaling) %>%
#     mutate(proapop_diff = (proliferation_plos - apoptosis_msigdbhallmark) / scaling) %>%
#     mutate(winapop_diff = (win - apoptosis_msigdbhallmark) / scaling) %>% 
#     mutate(winloseapop_diff = (win - apoptosis_msigdbhallmark - lose) / scaling)

# # value_names <- colnames(mosta)[c(6, 12:ncol(mosta))]
# # t_mosta <- cbind(
# #     mosta[, !colnames(mosta) %in% value_names],
# #     exp(mosta[, colnames(mosta) %in% value_names])
# # )

# #################################################################################
# ### Cor of Marker scores ########################################################
# #################################################################################

# # seltime <- "E9.5"
# # seltissue <- "Brain"

# lapply(levels(mosta$timepoint), function(seltime) {  
#     sel_mosta <- mosta %>%
#         filter(timepoint == seltime)
    
#     # all tisses
#     sel_mosta %>%
#         select(-cell_name, -newy, -newx, -x, -y, -breaks_ratio, -annotation, -timepoint, -section, -time) %>% 
#         cor %>% 
#         # ggcorrplot::ggcorrplot(hc.order = T, type = "lower", lab = TRUE) + 
#         ggcorrplot::ggcorrplot(lab = TRUE) +
#         ggtitle(str_glue("{seltime}"))
#     ggsave(str_glue("fig/cor/cor_{seltime}.png"), width = 30, height = 40, units = "in", bg = "white")

#     # separate tissue
#     pbmcapply::pbmclapply(unique(sel_mosta$annotation), function(seltissue) {
#         #cat(str_glue("{seltime}_{seltissue}"), "\n")
#         sel_mosta %>%
#             filter(annotation == seltissue) %>% 
#             select(-cell_name, -newy, -newx, -x, -y, -breaks_ratio, -annotation, -timepoint, -section, -time) %>% 
#             cor %>% 
#             # ggcorrplot::ggcorrplot(hc.order = T, type = "lower", lab = TRUE) + 
#             ggcorrplot::ggcorrplot(lab = TRUE) +
#             ggtitle(str_glue("{seltime}_{seltissue}"))
#         ggsave(str_glue("fig/cor/cor_{seltime}_{seltissue}.png"), width = 30, height = 40, units = "in", bg = "white")
#     }, mc.cores = 10)
# })


#################################################################################
### summary of emtropy resutls ##################################################
#################################################################################

rc <- 2
# th <- 10; numbin <- 10; rc <- 2;
lapply(c(5, 10, 20), function(th) {
    lapply(c(5, 10), function(numbin) {
        cat(str_glue("rc{rc}_th{th}_numbin{numbin}"), "\n")

        selfiles <- list.files("res/res_entropy/", full.names = TRUE, pattern = str_glue("{rc}_{th}_{numbin}.rds"))

        res <- pbmcapply::pbmclapply(selfiles, function(sel) {
            readRDS(sel)
        }, mc.cores = 20) %>% rbindlist()
        res$rc <- rc
        res$th <- th
        res$numbin <- numbin

        saveRDS(res, str_glue("res/summary_entropy_{rc}_{th}_{numbin}.rds"))
    })
})