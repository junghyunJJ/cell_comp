rm(list = ls())
library(tidyverse)
library(data.table)
library(tidytext)

library(jjutil)
library(openxlsx)

library(fgsea)

bodp_v1 <- fread("res/res_DAVID_mouseHALLMARK_MYC_TARGETS/GOBP_mouseHALLMARK_MYC_TARGETS_V1.txt")
sig_bodp_v1 <- bodp_v1 %>%
    filter(FDR < 0.05) %>%
    arrange(PValue) %>%
    mutate(Term = sub("(GO:.*)~(.*)", "\\2 (\\1)", Term))


#  (194 genes)
g1 <- ggplot(sig_bodp_v1[1:4, ], aes(reorder(Term, -FDR), -log10(FDR))) +
    geom_bar(stat = "identity", color = "black", fill = "firebrick3") +
    coord_flip() +
    theme_bw() +
    xlab("") +
    ylab("-log10(pvalue)") +
    ggtitle("MYC targets V1\n(HALLMARK)") +
    geom_hline(yintercept = -log10(0.01), col = "black", linetype = "dashed") +
    theme(
        axis.text = element_text(color = "black"),
        text = element_text(color = "black")
    )


bodp_v2 <- fread("res/res_DAVID_mouseHALLMARK_MYC_TARGETS/GOBP_mouseHALLMARK_MYC_TARGETS_V2.txt")
sig_bodp_v2 <- bodp_v2 %>%
    filter(FDR < 0.05) %>%
    arrange(PValue) %>%
    mutate(Term = sub("(GO:.*)~(.*)", "\\2 (\\1)", Term))

#  (58 genes)
g2 <- ggplot(sig_bodp_v2[1:4, ], aes(reorder(Term, -FDR), -log10(FDR))) +
    geom_bar(stat = "identity", color = "black", fill = "darkgreen") +
    coord_flip() +
    theme_bw() +
    xlab("") +
    ylab("-log10(pvalue)") +
    ggtitle("MYC targets V2\n(HALLMARK)") +
    geom_hline(yintercept = -log10(0.01), col = "black", linetype = "dashed") +
    theme(
        axis.text = element_text(color = "black"),
        text = element_text(color = "black")
    )

cowplot::plot_grid(g1, g2)
ggsave("fig/GOBP_MYC_TARGETS.png", height = 3, width = 10)
