rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)

library(jjutil)
library(openxlsx)

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

inner_cal_entropy_freq <- function(freq) {
    h <- -sum(ifelse(freq > 0, freq * log2(freq), 0))
    return(h)
}

cal_entropy_freq <- function(exp) {
    exp <- exp[apply(is.na(exp), 1, sum) == 0, ]
    probs <- t(t(exp) / apply(exp, 2, sum))
    entropy <- apply(probs, 2, inner_cal_entropy_freq)
    return(entropy)
}

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


# 2. analysis
dat <- NormalizeData(dat, verbose = FALSE)
dat <- FindVariableFeatures(dat, verbose = FALSE)
dat <- ScaleData(dat, features = rownames(dat), verbose = FALSE)
dat <- RunPCA(dat, verbose = FALSE)
dat <- FindNeighbors(dat, verbose = FALSE)
dat <- FindClusters(dat, verbose = FALSE)
dat


# 3. cal gene score
# 3-1. cal gene score using comp gene and hallmark

geneset <- readRDS("data/genesets.rds")
for (i in seq_len(length(geneset))) {
    dat <- AddModuleScore(
        dat,
        features = geneset[i], 
        ctrl = length(geneset[i]), 
        name = names(geneset[i])
    ) 
}
colnames(dat@meta.data) <- sub("(.*)1$", "\\1", colnames(dat@meta.data))

# diff values 
scaling <- 1
dat@meta.data <- dat@meta.data %>%
    mutate(winlose_diff = (win - lose) / scaling) %>%
    mutate(proapop_diff = (proliferation_plos - apoptosis_msigdbhallmark) / scaling) %>%
    mutate(winapop_diff = (win - apoptosis_msigdbhallmark) / scaling) %>% 
    mutate(winloseapop_diff = (win - apoptosis_msigdbhallmark - lose) / scaling)

# ratio value from Patric
res_cal_ratio <- cal_ratio(dat@assays$RNA@data, geneset$win, geneset$lose)
dat@meta.data <- dat@meta.data %>%
    mutate(ratio = res_cal_ratio$ratio) %>%
    mutate(breaks_ratio = res_cal_ratio$breaks)


# 4. cal entropy
# we focued on the normalized data
exp <- GetAssayData(dat)

# cal entropy
i <- 1
for (i in seq_len(length(geneset))) {
    print(i)
    sel_exp <- as.data.frame(exp)[geneset[[i]], ]
    entropy <- cal_entropy_freq(sel_exp)
    dat@meta.data[str_glue(names(geneset)[i], "_entropy")] <- entropy
}

# diff values 
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