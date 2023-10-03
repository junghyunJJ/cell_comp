
import warnings

import pandas as pd
import scanpy as sc
import moscot as mt

import re

warnings.simplefilter("ignore", UserWarning)


def grep(pattern, list):
    matched_lines = [line for line in list if re.search(pattern, line)]
    return matched_lines


allfiles = list(pd.read_table("allfiles_h5ad", header=None)[0])
allfiles = grep("MOSTA", allfiles)

for selfile in allfiles:
    print("########################")
    print(f"data/fromweb/{selfile}")
    print("########################")
    
    rawadata = sc.read(f"data/fromweb/{selfile}")
    rawadata.X = rawadata.layers["count"].copy()
    adata = rawadata.copy()

    # basic filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Remove cells that have too many mitochondrial genes expressed
    adata.var["mt"] = adata.var_names.str.startswith("mt-")  # annotate the group of mitochondrial genes as "mt" (#13 genes.. )
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    #sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
    adata = adata[adata.obs.pct_counts_mt < 5, :]

    # Normalization
    sc.pp.normalize_total(adata)

    # Logarithmize the data
    sc.pp.log1p(adata)

    # Identify highly-variable genes
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, subset=False)

    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # cal gene score
    # 1. comp gene
    win_marker=list(pd.read_csv("data/win_mouse.txt", header=None)[0])
    sc.tl.score_genes(adata,  win_marker, score_name="win", ctrl_size=len(win_marker))

    lose_marker=list(pd.read_csv("data/lose_mouse.txt", header=None)[0])
    sc.tl.score_genes(adata,  lose_marker, score_name="lose", ctrl_size=len(lose_marker))


    # 2. proliferation and apoptosis
    proliferation_moscot = mt.utils.data.proliferation_markers("mouse")
    # proliferation_moscot=list(pd.read_csv("data/proliferation_mouse_moscot.txt", header=None)[0])
    sc.tl.score_genes(adata, proliferation_moscot, score_name="proliferation_moscot", ctrl_size=len(proliferation_moscot))

    apoptosis_moscot = mt.utils.data.apoptosis_markers("mouse")
    # apoptosis_moscot=list(pd.read_csv("data/apoptosis_mouse_moscot.txt", header=None)[0])
    sc.tl.score_genes(adata, apoptosis_moscot, score_name="apoptosis_moscot", ctrl_size=len(apoptosis_moscot))

    proliferation_plos=list(pd.read_csv("data/proliferation_mouse_plos.txt", header=None)[0])
    sc.tl.score_genes(adata, proliferation_plos, score_name="proliferation_plos", ctrl_size=len(proliferation_plos))

    apoptosis_msigdbhallmark=list(pd.read_csv("data/apoptosis_mouse_msigdbhallmark.txt", header=None)[0])
    sc.tl.score_genes(adata, apoptosis_msigdbhallmark, score_name="apoptosis_msigdbhallmark", ctrl_size=len(apoptosis_msigdbhallmark))


    # 3. essentail gene
    essential_all=list(pd.read_csv("data/essential_mouse_all.txt", header=None)[0])
    sc.tl.score_genes(adata, essential_all, score_name="essential_all", ctrl_size=len(essential_all))

    essential_both=list(pd.read_csv("data/essential_mouse_both.txt", header=None)[0])
    sc.tl.score_genes(adata, essential_both, score_name="essential_both", ctrl_size=len(essential_both))


    # 4. haploinsufficiency gene
    haploinsufficiency=list(pd.read_csv("data/haploinsufficiency_mouse.txt", header=None)[0])
    sc.tl.score_genes(adata, haploinsufficiency, score_name="haploinsufficiency", ctrl_size=len(haploinsufficiency))
    save_selfile = re.sub(".h5ad", "_renormalized.h5ad", selfile)    
    #import pdb;pdb.set_trace()

    adata.write_h5ad("data/renormalized/" + save_selfile)
