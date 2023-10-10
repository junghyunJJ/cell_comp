rm(list = ls())
library(tidyverse)
library(data.table)

library(jjutil)
library(optparse)



allfiles_meta <- list.files("data/renormalized", pattern = "csv", full.names = TRUE)
allfiles <- list.files("data/renormalized", pattern = "rds", full.names = TRUE)

genesets <- readRDS("data/genesets.rds")

sel <- 49
z <- 1


sel_geneset <- genesets[[z]]
seurat <- readRDS(allfiles[sel])

dat <- GetAssayData(seurat)
# dat <- seurat@assays$RNA@data
dat <- seurat@assays$RNA@counts

sel_dat <- as.data.frame(dat)[sel_geneset, ]
sel_dat <- sel_dat[apply(is.na(sel_dat), 1, sum) == 0, ]
sel_dat %>% dim
sel_dat %>% h
sel_dat %>% colSums()

entropyJJ <- function(freqs) {
    H <- -sum(ifelse(freqs > 0, freqs * log(freqs), 0))
    H <- H / log(2)
    return(H)
}
entropyJJ2 <- function(ptab) {
    H <- -sum(ifelse(ptab > 0, ptab * log2(ptab), 0))
    return(H)
}
    

probs <- t(t(sel_dat) / apply(sel_dat, 2, sum))
probs %>% colSums()
apply(probs, 2, entropyJJ) %>% summary
apply(probs, 2, entropyJJ2) %>% summary

apply(probs, 2, entropy::entropy.plugin, unit = "log2") %>% summary

# entropy <- -apply(probs * log2(probs) / log2(nrow(sel_dat)), 2, sum)
entropy <- -apply(probs * log2(probs) / log(2), 2, sum)
entropy %>% summary
# Generated from function body. Editing this file has no effect.

freqs <- probs[, 1]
freqs <- rep(1 / 6, 6)

entropyJJ(rep(1 / 6, 6))
entropyJJ2(rep(1 / 6, 6))
entropy::entropy.plugin(rep(1 / 6, 6), unit = "log2")
