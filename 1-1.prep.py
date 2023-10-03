
import warnings

import pandas as pd
import scanpy as sc

import re

import gseapy as gp

def grep(pattern, list):
    matched_lines = [line for line in list if re.search(pattern, line)]
    return matched_lines

warnings.simplefilter("ignore", UserWarning)

# read mouse hallmark
mh = gp.read_gmt(path="/home/jungj2/tools/msigdb/msigdb_v2023.1.Mm_files_to_download_locally/msigdb_v2023.1.Mm_GMTs/mh.all.v2023.1.Mm.symbols.gmt")

allfiles = list(pd.read_table("allfiles_h5ad", header=None)[0])
allfiles = grep("MOSTA", allfiles)

for rawselfile in allfiles:
    selfile = re.sub(".h5ad", "_renormalized.h5ad", rawselfile)    
    print(selfile)
    
    # adata = sc.read("data/renormalized/E15.5_E2S1.MOSTA_renormalized.h5ad")
    adata = sc.read("data/renormalized/" + selfile)
    adata.obs = adata.obs[['win', 'lose', 'proliferation_moscot', 'apoptosis_moscot', 'proliferation_plos', 'apoptosis_msigdbhallmark', 'essential_all', 'essential_both', 'haploinsufficiency']]

    for i in range(len(mh)):
        key, value = list(mh.items())[i]
        print(f"{i}: {key} ({len(value)})")
        
        sc.tl.score_genes(adata,  value, score_name=key, ctrl_size=len(value))
        selfile = re.sub(".h5ad", "_renormalized.h5ad", rawselfile)    
        save_selfile = re.sub(".h5ad", "_obs.csv", selfile)    
    
    #import pdb;pdb.set_trace()
    adata.obs['file'] = selfile
    adata.obs.to_csv("data/renormalized/" + save_selfile)
