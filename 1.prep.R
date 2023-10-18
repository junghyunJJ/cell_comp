rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)

library(jjutil)
library(openxlsx)

library(Seurat)
source("functions.R")

sel <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# sel <- 1

#############################################################################
#############################################################################
#############################################################################

allfiles <- readLines("allfiles_h5ad")
allfiles <- grep("MOSTA", allfiles, value = TRUE)
selfile <- allfiles[sel]

cat("#################################", sep = "\n")
cat(str_glue("data/fromweb/{selfile}"), sep = "\n")
cat("#################################", sep = "\n")

selname <- sub("(E.*)\\.MOSTA\\.h5ad", "\\1", selfile)
rawdat <- anndata::read_h5ad(str_glue("data/fromweb/{selfile}"))


# 1. data filtering 
# 1-1. create Seurat object
dat <- CreateSeuratObject(counts = t(as.matrix(rawdat$layers["count"])), project = selname)


# 1-2. add meta information
if (!all.equal(rownames(rawdat$obs), rownames(dat@meta.data))) {
    stop("PLEASE check data")
}

tmp <- unlist(strsplit(selname, "_"))
timepoint <- tmp[1]
time <- sub("E", "", timepoint) %>% as.numeric()
section <- sub(".MOSTA", "", tmp[2])

dat@meta.data <- dat@meta.data %>%
    mutate(annotation = rawdat$obs$annotation) %>%
    mutate(timepoint = timepoint) %>% 
    mutate(time = time) %>% 
    mutate(section = section)


# 1-3. add coord data
coord <- rawdat$obsm$spatial %>% as.data.frame %>% 
    set_colnames(c("x", "y"))
dat@meta.data <- cbind(dat@meta.data, coord)


# 1-4. basic filtering
dat <- CreateSeuratObject(GetAssayData(dat), min.cells = 3, min.features = 200, meta.data = dat@meta.data)
dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-")
dat <- subset(dat, subset = nCount_RNA > 200 & percent.mt < 5)
dat


# # 1-5. check the data with results from Scnapy
# scanpy <- anndata::read_h5ad(str_glue("data/oldrenormalized/{timepoint}_{section}.MOSTA_renormalized.h5ad"))
# if (nrow(scanpy$var) != length(rownames(dat))) {
#     stop("PLEASE check min.cells options for basic filtering")
# }

# if (nrow(scanpy$obs) != length(colnames(dat))) {
#     stop("PLEASE check min.features options or mt gene filtering for basic filtering")
# }


# 2. analysis ( we ONLY need the NormalizeData step at this point)
dat <- NormalizeData(dat, verbose = FALSE)
dat <- FindVariableFeatures(dat, verbose = FALSE)
dat <- ScaleData(dat, features = rownames(dat), verbose = FALSE)
dat <- RunPCA(dat, verbose = FALSE)
dat <- FindNeighbors(dat, verbose = FALSE)
dat <- FindClusters(dat, verbose = FALSE)
dat


# 3. cal gene score and entropy
# 3-1. cal gene score using comp gene and hallmark

geneset <- readRDS("data/genesets.rds")[1:2]
for (i in seq_len(length(geneset))) {
    dat <- AddModuleScore(
        dat,
        features = geneset[i], 
        ctrl = length(geneset[[i]]) / 3, 
        name = paste0("ori_", names(geneset[i]))
    ) 
}
colnames(dat@meta.data) <- sub("(.*)1$", "\\1", colnames(dat@meta.data))

# 3-2. cal entropy
for (i in seq_len(length(geneset))) {

    dat <- oldcal_genescore(
        seurat_obj = dat,
        geneset = geneset[[i]], 
        background = round(length(geneset[[i]]) / 3),
        name = paste0("old_", names(geneset[i]))
    ) 
    dat <- cal_genescore(
        seurat_obj = dat,
        geneset = geneset[[i]], 
        background = round(length(geneset[[i]])),
        name = names(geneset[i])
    ) 
    dat <- cal_entropy(
        seurat_obj = dat,
        geneset = geneset[[i]],
        background = round(length(geneset[[i]])),
        name = names(geneset[i])
    )
    dat <- cal_entropy(
        seurat_obj = dat,
        geneset = geneset[[i]], 
        background = round(length(geneset[[i]]) / 3),
        name = paste0(names(geneset[i]), 2)
    ) 
}
# dat@meta.data %>% head
# cor(dat@meta.data$ori_win, dat@meta.data$old_win_genescore)
# cor(dat@meta.data$ori_win, dat@meta.data$win_genescore)
# cor(dat@meta.data$win_genescore, dat@meta.data$win_entropy, use = "complete")
# plot(dat@meta.data$win_genescore, dat@meta.data$win_entropy)

# cor(dat@meta.data$win_entropy, dat@meta.data$win2_entropy, use = "complete")
# plot(dat@meta.data$win_entropy, dat@meta.data$win2_entropy)


# genescore diff values 
scaling <- 1
dat@meta.data <- dat@meta.data %>%
    mutate(winlose_diff = (win_genescore - lose_genescore) / scaling) %>%
    mutate(proapop_diff = (proliferation_plos_genescore - apoptosis_msigdbhallmark_genescore) / scaling) %>%
    mutate(winapop_diff = (win_genescore - apoptosis_msigdbhallmark_genescore) / scaling) %>% 
    mutate(winloseapop_diff = (win_genescore - apoptosis_msigdbhallmark_genescore - lose_genescore) / scaling)

# ratio value from Patric
res_cal_ratio <- cal_ratio(dat@assays$RNA@data, geneset$win, geneset$lose)
dat@meta.data <- dat@meta.data %>%
    mutate(ratio = res_cal_ratio$ratio) %>%
    mutate(breaks_ratio = res_cal_ratio$breaks)


# entropy diff values 
dat@meta.data <- dat@meta.data %>%
    mutate(winlose_diff_entropy = (win_entropy - lose_entropy) / scaling) %>%
    mutate(proapop_diff_entropy = (proliferation_plos_entropy - apoptosis_msigdbhallmark_entropy) / scaling) %>%
    mutate(winapop_diff_entropy = (win_entropy - apoptosis_msigdbhallmark_entropy) / scaling) %>% 
    mutate(winloseapop_diff_entropy = (win_entropy - apoptosis_msigdbhallmark_entropy - lose_entropy) / scaling)


# 5. save file
dat@meta.data <- dat@meta.data %>%
    rownames_to_column("cell_name") %>%
    mutate(timepoint = timepoint) %>%
    mutate(time = time) %>%
    mutate(section = section) %>%
    as_tibble()

write_csv(dat@meta.data, str_glue("data/renormalized/{timepoint}_{section}.MOSTA_renormalized_obs.csv"))
saveRDS(dat, str_glue("data/renormalized/{timepoint}_{section}.MOSTA_renormalized.rds"))