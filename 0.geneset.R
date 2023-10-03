rm(list = ls())
library(tidyverse)
library(data.table)

library(jjutil)
library(openxlsx)

convert_mouse_to_human <- function(gene_list) { 
    output <- c()
    mouse_human_genes <- read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

    for (gene in gene_list) {
        class_key <- (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "mouse, laboratory"))[['DB.Class.Key']]
        if( !identical(class_key, integer(0)) ) {
            human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
            for(human_gene in human_genes) {
            output = rbind(c(gene, human_gene), output)
            }
        }
    }
    return (output)
}

convert_human_to_mouse <- function(gene_list) {
  output = c()
  mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  for(gene in gene_list) {
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "human"))[['DB.Class.Key']]
    if( !identical(class_key, integer(0)) ) {
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      for(human_gene in human_genes) {
        output = rbind(c(gene, human_gene), output)
      }
    }
  }
  return (output)
}


#################################################################
### Comp gene ################################################### 
#################################################################

rawcompgene <- openxlsx::read.xlsx("data/comp_gene.xlsx")
win_mouse <- rawcompgene %>% filter(type == "win") %>% select(mouse_sym) %>% pull %>% sort
lose_mouse <- rawcompgene %>% filter(type == "lose") %>% select(mouse_sym) %>% pull %>% sort
writeLines(win_mouse, "data/win_mouse.txt")
writeLines(lose_mouse, "data/lose_mouse.txt")


#################################################################
### Apoptosis ################################################### 
#################################################################

# The gene list for apoptosis is available from msigdb.
# https://www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/HALLMARK_APOPTOSIS.html

mh <- read.gmt("/home/jungj2/tools/msigdb/msigdb_v2023.1.Mm_files_to_download_locally/msigdb_v2023.1.Mm_GMTs/mh.all.v2023.1.Mm.symbols.gmt")
apoptosis_mouse <- mh$HALLMARK_APOPTOSIS %>% sort
writeLines(apoptosis_mouse, "data/apoptosis_mouse_msigdbhallmark.txt")
apoptosis_mouse %>% ll # 161

apoptosis_moscot_mouse <- readLines("data/apoptosis_mouse_moscot.txt")
apoptosis_moscot_mouse %>% ll # 193

# NOTE!! only small fraction fo genes are overlapped between msigdb and moscot list because the moscot list come from Hallmark P53 Pathway,msigdb.
intersect(apoptosis_mouse, apoptosis_moscot_mouse) %>% ll # 21


# h <- read.gmt("/home/jungj2/tools/msigdb/msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/h.all.v2023.1.Hs.symbols.gmt")
# apoptosis_human <- h$HALLMARK_APOPTOSIS %>% sort
# apoptosis_human %>% ll # 161
# writeLines(apoptosis_human, "data/proliferation_human.txt")

#################################################################
### Proliferation ###############################################
#################################################################

# 1. The gene list for apoptosis is available from GO.
# cell population proliferation (GO:0008283) / https://amigo.geneontology.org/amigo/term/GO:0008283
# https://www.informatics.jax.org/go/term/GO:0008283

# NOTE! we onlt focused on the Experimental evidence codes
# https://geneontology.org/docs/guide-go-evidence-codes/

sel_gocode <- c("EXP","IDA","IPI","IMP","IGI","IEP","HTP","HDA","HMP","HGI","HEP")

# Proliferation
all_GO0008283_mouse_go <- read.xlsx("data/GO_term_summary_20230907_023720.xlsx") %>% as_tibble()
all_GO0008283_mouse_go %>% dim # 10426

f_GO0008283_mouse_go <- all_GO0008283_mouse_go %>% filter(Evidence %in% sel_gocode)
f_GO0008283_mouse_go %>% dim # 5862
f_GO0008283_mouse_go$Evidence %>% table
proliferation_mouse_go <- f_GO0008283_mouse_go$Symbol %>% unique %>% sort
proliferation_mouse_go

# 2. msigdb
all_msigdb <- jjutil::read.gmt("/home/jungj2/tools/msigdb/msigdb_v2023.1.Mm_files_to_download_locally/msigdb_v2023.1.Mm_GMTs/m5.go.bp.v2023.1.Mm.symbols.gmt")
tmp_msigdb <- grep("GOBP_.*_PROLIFERATION$",names(all_msigdb), value = T)
proliferation_msigdb <- tmp_msigdb[-grep("GOBP_NEGATIVE_REGULATION_OF_|GOBP_POSITIVE_REGULATION_OF_|GOBP_REGULATION_OF_", tmp_msigdb)] %>% sort

proliferation_mouse_msigdb <- all_msigdb[names(all_msigdb) %in% proliferation_msigdb] %>% unlist %>% unique %>% sort
proliferation_mouse_msigdb %>% ll # 1456

# 3. moscot
proliferation_mouse_moscot <- readLines("data/proliferation_mouse_moscot.txt")
proliferation_mouse_moscot %>% ll # 98


# 4. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010604#sec018
proliferation_plos <- read.xlsx("data/proliferation_gene.xlsx")
proliferation_mouse_plos <- proliferation_plos$mouse_sym
proliferation_mouse_plos %>% ll # 25
writeLines(proliferation_mouse_plos, "data/proliferation_mouse_plos.txt")


#################################################################
### Haploinsufficiency ########################################## 
#################################################################
# https://search.clinicalgenome.org/kb/downloads

raw_clinGen <- fread("data/Clingen-Curation-Activity-Summary-Report-2023-09-07.csv",skip = 3)
raw_clinGen$haploinsufficiency_score <- sub("(.*) \\([0-9]+/[0-9]+/[0-9]+\\)$", "\\1",raw_clinGen$dosage_haploinsufficiency_assertion)
raw_clinGen$haploinsufficiency_score %>% table %>% data.frame

haploinsufficiency_gene <- raw_clinGen %>% 
    filter(haploinsufficiency_score == "3 - Sufficient Evidence for Haploinsufficiency") %>% 
    select(gene_symbol) %>% pull %>% unique %>% sort
haploinsufficiency_gene %>% ll

haploinsufficiency_mouse <- convert_human_to_mouse(haploinsufficiency_gene)[,2] %>% unique %>% sort
writeLines(haploinsufficiency_mouse, "data/haploinsufficiency_mouse.txt")



# seuratObject <- LoadH5Seurat("data/fromweb/E11.5_E1S3.MOSTA.h5seurat")
# seuratObject  %>% str
# seuratObject@reductions$spatial@cell.embeddings %>% dim
# seuratObject@meta.data %>% dim

# dat <- data.frame(seuratObject@reductions$spatial@cell.embeddings, annotation = seuratObject@meta.data$annotation) %>% as_tibble
# dat

# #col <- WGCNA::labels2colors(compgene_mosta$annotation) %>% unique

# ggplot(dat, aes(spatial_1, spatial_2, col = annotation)) +
#   geom_point(size = 0.5) +
#   # theme_void()+
#   theme_bw()