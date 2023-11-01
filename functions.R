library(tidyverse)
library(data.table)
library(magrittr)

library(jjutil)

oldcal_genescore <- function(seurat_obj, geneset, background = NULL, nbin = 24, name = NULL, seed = 1, rawgenescore = FALSE) {
    set.seed(seed)

    # emparically the length(geneset) / 3 looks resonable
    if (is.null(background)) {
        background <- length(geneset)
    }

    if (is.null(name)) {
        name <- "geneset"
    }

    # we use the expression data based on the seurat_obj
    exp <- GetAssayData(seurat_obj)

    # filter geneset based on expression data
    geneset <- geneset[geneset %in% rownames(exp)]
    
    # All analyzed features are binned based on averaged expression
    dat_avg <- exp %>%
        rowMeans() %>%
        sort()
    
    # cut data
    dat_cut <- cut_number(dat_avg + rnorm(n = length(dat_avg)) / 1e+30, n = nbin, labels = FALSE, right = FALSE)
    names(dat_cut) <- names(dat_avg)

    # the background genes are randomly selected from each bin.
    background_genes <- lapply(seq_along(geneset), function(j) {
        sample(dat_cut[which(dat_cut == dat_cut[geneset[j]])], size = background, replace = FALSE) %>% names
    }) %>% unlist %>% unique

    # cal score
    background_scores <- Matrix::colMeans(exp[background_genes, ])
    scores <- Matrix::colMeans(exp[geneset, ])
    
    seurat_obj@meta.data[paste0(name, "_genescore")] <- scores - background_scores

    if (rawgenescore) {
        seurat_obj@meta.data[paste0(name, "_rawgenescore")] <- scores
    }
    
    return(seurat_obj)
}

cal_genescore <- function(seurat_obj, geneset, background = NULL, nbin = 24, name = NULL, seed = 1, rawgenescore = FALSE) {
    set.seed(seed)

    # emparically the length(geneset) / 3 looks resonable
    if (is.null(background)) {
        background <- length(geneset)
    }

    if (is.null(name)) {
        name <- "geneset"
    }

    # we use the expression data based on the seurat_obj
    exp <- GetAssayData(seurat_obj) %>% as.matrix()

    # filter geneset based on expression data
    geneset <- geneset[geneset %in% rownames(exp)]

    # All analyzed features are binned based on averaged expression
    dat_avg <- apply(exp, 1, median) %>%
        sort()

    # cut data
    dat_cut <- cut_number(dat_avg + rnorm(n = length(dat_avg)) / 1e+30, n = nbin, labels = FALSE, right = FALSE)
    names(dat_cut) <- names(dat_avg)

    # the background genes are randomly selected from each bin.
    background_genes <- lapply(seq_along(geneset), function(j) {
        sample(dat_cut[which(dat_cut == dat_cut[geneset[j]])], size = background, replace = FALSE) %>% names()
    })

    # Calculate the background_genescore using median value of each selected bin
    tmp_background_scores <- lapply(background_genes, function(sel) {
        # data.frame(Matrix::colMeans(exp[sel, ])) %>% t %>% as.data.frame()
        data.frame(apply(exp[sel, ], 2, median)) %>%
            t() %>%
            as.data.frame()
    }) %>% rbindlist()
    background_scores <- apply(tmp_background_scores, 2, median, na.rm = TRUE)

    # cal score
    # background_scores <- Matrix::colMeans(exp[unique(unlist(background_genes)), ])
    scores <- Matrix::colMeans(exp[geneset, ])

    seurat_obj@meta.data[paste0(name, "_genescore")] <- scores - background_scores

    if (rawgenescore) {
        seurat_obj@meta.data[paste0(name, "_rawgenescore")] <- scores
    }

    return(seurat_obj)
}

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
    exp <- as.matrix(exp)
    exp <- exp[apply(is.na(exp), 1, sum) == 0, ]
    probs <- t(t(exp) / apply(exp, 2, sum))
    entropy <- apply(probs, 2, inner_cal_entropy_freq)
    return(entropy)
}

cal_gene_entropy_freq <- function(exp) {
    exp <- as.matrix(exp)
    exp <- exp[apply(is.na(exp), 1, sum) == 0, ]
    probs <- exp / apply(exp, 1, sum)
    entropy <- apply(probs, 1, inner_cal_entropy_freq)
    return(entropy)
}

cal_entropy <- function(seurat_obj, geneset, background = NULL, nbin = 24, name = NULL, seed = 1, rawgentropy = FALSE) {
    
    set.seed(seed)

    # emparically the length(geneset) / 3 looks resonable
    if (is.null(background)) {
        background <- length(geneset)
    }

    if (is.null(name)) {
        name <- "geneset"
    }

    # we use the expression data based on the seurat_obj
    exp <- GetAssayData(seurat_obj)

    # filter geneset based on expression data 
    geneset <- geneset[geneset %in% rownames(exp)]

    # All analyzed features are binned based on gene based entropy
    gene_entropy <- cal_gene_entropy_freq(exp) %>% sort
    
    # the background genes are randomly selected from each bin.
    dat_cut <- cut_number(gene_entropy + rnorm(n = length(gene_entropy)) / 1e+30, n = nbin, labels = FALSE, right = FALSE)
    names(dat_cut) <- names(gene_entropy)

    # randomly select background genes
    background_genes <- lapply(seq_along(geneset), function(j) {
        sample(dat_cut[which(dat_cut == dat_cut[geneset[j]])], size = background, replace = FALSE) %>% names
    })
    
    # Calculate the background_entropy using median value of each selected bin
    tmp_background_entropy <- lapply(background_genes, function(sel) {
        data.frame(cal_entropy_freq(exp[sel, , drop  = FALSE])) %>% t %>% as.data.frame()
    }) %>% rbindlist
    background_entropy <- apply(tmp_background_entropy, 2, median, na.rm = TRUE)

    entropy <- cal_entropy_freq(exp[geneset, ])

    seurat_obj@meta.data[paste0(name, "_entropy")] <- entropy - background_entropy

    if (rawgentropy) {
        seurat_obj@meta.data[paste0(name, "_rawentropy")] <- entropy
    }

    return(seurat_obj)
}

# dat "trans_dat" function is exactly same as "scaleTR" function in "cdfquantreg" r pakcage
trans_dat <- function(y, high = NULL, low = NULL, data = NULL, N = NULL, scale = 0.5) {
    if (!is.null(data)) {
        yname <- deparse(substitute(y))
        y <- data[, yname]
    }
    if (is.null(high)) {
        high <- max(y, na.rm = TRUE)
    }
    if (is.null(low)) {
        low <- min(y, na.rm = TRUE)
    }
    if (is.null(N)) {
        N <- length(na.omit(y))
    }
    y0 <- y[!is.na(y)]
    y1 <- (y0 - low) / (high - low)
    y2 <- (y1 * (N - 1) + scale) / N
    y3 <- y
    y3[!is.na(y)] <- y2
    if (is.null(data)) {
        return(y3)
    } else {
        names <- paste(yname, "old", sep = "")
        data[, names] <- y
        data[, yname] <- y3
        return(data)
    }
}

# Matrix::sparseMatrix(i = dat$x, j = dat$y, x = dat$value) %>% as.matrix
make_mat <- function(dat) {
    dat <- dat %>%
        dplyr::select("x", "y", "value", "annotation", "cellname") %>%
        mutate(cellname = sub("(.*)-.*-.*", "\\1", cellname))
    
    xlim <- round(max(dat$x))
    ylim <- round(max(dat$y))
    
    # make null matrix  
    mat_dat <- matrix(rep(NA, xlim * ylim), nrow = xlim, ncol = ylim)
    mat_anno <- matrix(rep(NA, xlim * ylim), nrow = xlim, ncol = ylim)
    mat_name <- matrix(rep(NA, xlim * ylim), nrow = xlim, ncol = ylim)

    # i<-1
    for (i in seq_len(nrow(dat))) {
        z <- dat[i, ]
        mat_dat[z[[1]], z[[2]]] <- z[[3]]
        mat_anno[z[[1]], z[[2]]] <- as.character(z[[4]])
        mat_name[z[[1]], z[[2]]] <- as.character(z[[5]])        
        # NOTE!!!!!!! filename is next step
        # mat_anno[z[[1]], z[[2]]] <- paste0(as.character(z[[4]]), "/", as.character(z[[5]])) 
    }
    return(list(mat_dat = mat_dat, anno_dat = mat_anno, name_dat = mat_name))
}

# make bin matrix (default <- row = 2; column = 2)
bin_mat <- function(M, r = 2, c = 2) {
    nr <- ceiling(nrow(M) / r)
    nc <- ceiling(ncol(M) / c)
    newM <- matrix(NA, nr * r, nc * c)
    newM[seq_len(nrow(M)), seq_len(ncol(M))] <- M
    
    div_k <- kronecker(matrix(seq_len(nr * nc), nr, byrow = TRUE), matrix(1, r, c))
    matlist <- split(newM, div_k)
    N <- length(matlist)
    mats <- unlist(matlist)
    dim(mats) <- c(r, c, N)
    res <- plyr::alply(mats, 3)
    return(res)
    return(res)
}

# split marix based on the specific bin (default <- row = 2; column = 2)
grid_mat <- function(matdat, r = 2, c = 2, complete = FALSE) {
    g_matdat <- bin_mat(matdat$mat_dat, r, c)
    anno_matdat <- bin_mat(matdat$anno_dat, r, c)
    name_matdat <- bin_mat(matdat$name_dat, r, c)
    
    if (complete) {
        check_na <- lapply(anno_matdat, function(z) {sum(is.na(z))}) == 0
        g_matdat <- g_matdat[check_na]
        anno_matdat <- anno_matdat[check_na]
        name_matdat <- name_matdat[check_na]
    }

    anno <- lapply(anno_matdat, function(sel_anno) {
        sel_anno <- sel_anno %>% 
            as.vector() %>% 
            unique()
      
        sel_anno <- sel_anno[!is.na(sel_anno)]
        paste0(sort(sel_anno), collapse = ";")
    }) %>% unlist

    names(g_matdat) <- anno
    names(name_matdat) <- anno
    return(list(matdat = g_matdat, cellname = name_matdat))
}

# make cell type com
# make_celltype_comb <- function(celltype, includehomo = FALSE, exclude = NULL) {
#     # remove specific celltype
#     if (!is.null(exclude)) {
#         celltype <- celltype[!celltype %in% exclude]
#     }

#     res_comb <- apply(expand.grid(celltype, celltype), 1, function(ss) {
#         ss %>%
#             as.vector() %>%
#             unique() %>%
#             sort() %>%
#             paste(collapse = ";")
#     }) %>% unique()

#     # if you want to include homo
#     if (!includehomo) {
#         res_comb <- grep(";", res_comb, value = TRUE)
#     }
#     return(res_comb)
# }

make_celltype_comb <- function(g_matdat, th = 10, includehomo = TRUE, exclude = NULL) {
    celltype_comb <- g_matdat$matdat %>%
        names() %>%
        table() %>%
        sort()
    
    # remove specific celltype
    if (!is.null(exclude)) {
        celltype_comb <- celltype_comb[!grepl(paste(exclude, collapse = "|"), names(celltype_comb))]        
    }

    celltype_comb <- celltype_comb[celltype_comb > th] %>% names
    celltype_comb <- celltype_comb[celltype_comb != ""]
    if (!includehomo) {
        celltype_comb <- grep(";", celltype_comb, value = TRUE)
    }
    return(celltype_comb)
}

summary_score <- function(g_matdat, celltype_comb, name = NULL) {
    res_celltype_comb <- lapply(celltype_comb, function(sel) {
        # print(sel)
        idx <- names(g_matdat$matdat) == sel
        raw_sel_dat <- g_matdat$matdat[idx]
        raw_sel_name <- g_matdat$cellname[idx]

        # data
        sel_dat <- raw_sel_dat %>%
            unlist() %>%
            as.numeric()
        
        # name
        sel_name <- raw_sel_name %>%
            unlist() %>%
            as.character()
        
        # summary
        res <- data.frame(
            cellname = sel_name, 
            celltype = sel,
            value = sel_dat, 
            # value = sel_dat[!is.na(sel_dat)], 
            n_grid = length(raw_sel_dat)
        )
        
        # remove NAs
        res <- res %>%
            filter(!is.na(cellname)) 

        # rename for value column
        colnames(res)[colnames(res) == "value"] <- paste(name)
        res
        
    }) %>% rbindlist()
    return(res_celltype_comb)
}

sel_surrounding <- function(homo_res, type = "apoptosis", score = "genescore", n_top = 10) {

    sel_type <- grep(type, colnames(homo_res), value = TRUE) %>%
        grep(score, ., value = TRUE)
    cat("type:", sel_type, " / # of center cell:", n_top, "\n")

    sel_dat <- homo_res %>%
        arrange(desc(get(sel_type))) %>%
        top_n(get(sel_type), n = n_top)

    res <- lapply(seq_len(nrow(sel_dat)), function(i) {
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
        s_dat <- homo_res %>%
            unite("coord", y:x, sep = "_", remove = FALSE) %>%
            filter(coord %in% s_coord) %>%
            select(-coord, -x, -y, -annotation, -n_grid)
        # browser()

        # check celltype: we only focued on the homo
        check_celltype <- c(c_dat$celltype, s_dat$celltype) %>% unique()
        
        # we skip res if there are NA in >50%
        if (length(check_celltype) == 1 & sum(!is.na(s_dat$winlose_diff_entropy)) > 3 & sum(!is.na(s_dat$winlose_diff_genescore)) > 3) {
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
        } else {
            summary_dat <- NA
        }
        summary_dat
    })
    return(res)
}


# # cal entropy value using discretize results
# cal_entropy_value <- function(dat, numbin = 10, unit = "log2") {
#     entropy::discretize(dat, numBins = numbin, r = c(0, 1)) %>%
#         entropy::entropy.empirical(unit = unit)
# }

# # cal entropy based on location (i.e., edge of the tissue)
# cal_entropy_whole_tissue <- function(dat, celltype, numbin = 10, unit = "log2") {
#     res_entropy <- lapply(celltype, function(sel) {
#         sel_values <- dat %>%
#             filter(annotation == sel) %>%
#             select(value) %>%
#             pull()
#         cal_entropy_value(sel_values, numbin = numbin, unit = unit)
#     }) %>% unlist
#     return(data.frame(celltype = celltype, entropy = res_entropy))
# }

# # cal entropy based on location (i.e., edge of the tissue)
# cal_entropy <- function(s_matdat, celltype_comb, numbin = 10, unit = "log2") {
#     res_entropy <- lapply(celltype_comb, function(sel_comb) {
#         sel_location <- s_matdat[names(s_matdat) == sel_comb] %>% unlist()      
#         sel_location <- sel_location[!is.na(sel_location)]
#         cal_entropy_value(sel_location, numbin = numbin, unit = unit)
#     }) %>% unlist
#     return(data.frame(celltype = celltype_comb, entropy = res_entropy))
# }