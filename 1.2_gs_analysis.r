# This script performs gene-set analyses based on differential gene expression result obtained from 1.1_DESeq2.r
# Output from this script corresponds to Figure 4B, C, Extended Figure 7F, I
library(dplyr)
library(clusterProfiler)
library(simplifyEnrichment)
library(GOSemSim)

# Read DESeq2 result============================================
res_compile = readRDS('./result/res_compile.RDS') #modify based on path of res_compile

# Gene ontology over-representation analysis================================
go_res = lapply(1:3, function(n) {
  
  res = res_compile[[n]]
  
  # Munge summary statistics
  res_ordered = res %>%
    as.data.frame() %>%
    arrange(padj) %>%
    rownames_to_column('ENTREZID')
  res_ordered = res_ordered[!is.na(res_ordered$pvalue) & !is.na(res_ordered$padj), ]
  up = subset(res_ordered, log2FoldChange > 0)
  down = subset(res_ordered, log2FoldChange < 0)
  up_sig = up[up$padj <= 0.05, ]
  down_sig = down[down$padj <= 0.05, ]
  
  
  sig_list = list(up = up_sig, down = down_sig)
  go_res = lapply(sig_list, function(i) {
    
    res = enrichGO(i$ENTREZID,
                   OrgDb = org.Mm.eg.db, #If human, use org.Hs.eg.db
                   keyType = 'ENTREZID',
                   ont = 'BP',
                   readable = F) 

  }) %>% setNames(names(sig_list))
  
  return(go_res)
  
}) %>% setNames(names(res_compile)[1:3])

# Gene set enrichment analysis================================================
destfile = "m5.all.v2022.1.Mm.entrez.gmt"
c5 = read.gmt(destfile)
c5 = c5[grep('GOBP', c5$term), ]
gsea_res = lapply(1:3, function(n) {
  
  res = res_compile[[n]]
  
  # Munge summary statistics
  res_ordered = res %>%
    as.data.frame() %>%
    arrange(padj) %>%
    rownames_to_column('ENTREZID')
  res_ordered = res_ordered[!is.na(res_ordered$pvalue) & !is.na(res_ordered$padj), ]
  res_ordered = res_ordered %>%
    mutate(rank = ifelse(log2FoldChange > 0, -log(pvalue), -log(pvalue)*-1))
  
  # Create rank input
  rnk = res_ordered$rank
  names(rnk) = res_ordered$ENTREZID
  rnk = sort(rnk, decreasing = TRUE)
  
  # GSEA
  c5_enrich = GSEA(rnk, TERM2GENE = c5, pvalueCutoff = 1)
  
}) %>% setNames(names(res_compile)[1:3])


# Semantic clustering===============================================
# Load requisite files
mmGO = godata('org.Mm.eg.db', ont = "BP")
go_id = read.delim('go_id.txt', header = F)
go_id$V2 = gsub(" ", "", go_id$V2) %>%
  tolower()

# Munge
c5_res = gsea_res[[1]][, c(5,7)] %>%
  rownames_to_column('ID')
c5_id = lapply(c('up', 'down'), function(i) {
  
  if (i == 'up') c5_id = subset(c5_res, NES > 0) 
  if (i == 'down') c5_id = subset(c5_res, NES < 0) 
  c5_id = c5_id %>%
    extract(1) %>%
    as.data.frame() %>%
    set_colnames("ID")
  c5_id$ID = gsub("GOBP_|_", "", c5_id$ID) %>%
    tolower()
  c5_id = left_join(c5_id, go_id, by = c("ID" = "V2")) 
  c5_id = na.omit(c5_id$V1)
  
}) %>% setNames(c('up', 'down'))

# Cluster GO terms
up_gosim = GO_similarity(c5_id$up, ont = 'BP')
dn_gosim = GO_similarity(c5_id$down, ont = 'BP')
simplify_up = simplifyGO(up_gosim, plot = TRUE)
simplify_dn = simplifyGO(dn_gosim, plot = TRUE)
