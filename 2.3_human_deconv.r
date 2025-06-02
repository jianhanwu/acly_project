# This script performs differential gene expression analysis between disease states using CEL files obtained from GSE164760
# Output from this script corresponds to Figure 5N, Extended Figure 12A, B, C
library(dplyr)
library(affy)
library(magrittr)
library(tidyverse)
BiocManager::install('metapredict', site_repository = 'https://hugheylab.github.io/drat/')
library(metapredict)
installCustomCdfPackages("hgu219hsentrezgcdf", ver = 25)
library(hgu219hsentrezgcdf)
library(immunedeconv)

# Load gene ID
id = read.delim('hgu219_entrez_symbol.txt', header = T)

# Get file path for RMA data
dir_path = '/path/to/GSE164760/'

# Read
GSE164760 = ReadAffy(celfile.path = dir_path, cdfname="hgu219hsentrezg")

# Create expression file
norm.GSE164760 = rma(GSE164760)
GSE164760.exp = as.data.frame(exprs(norm.GSE164760))
rownames(GSE164760.exp) = gsub('_at', '', rownames(GSE164760.exp))

# Exponentiate logarithmic scale expression
GSE164760.norm = 2^GSE164760.exp %>%
  rownames_to_column('ENTREZID') %>%
  cbind(Symbol = id[match(.$ENTREZID, id$ENTREZID), ]$Symbol) %>%
  filter(!is.na(Symbol)) %>%
  column_to_rownames('Symbol')
GSE164760.norm = GSE164760.norm[, -1]

# Correlation between ACLY and B cell markers====================================================
gs_gepliver = read.csv('./gepliver_cell_annot.csv')
gs_list_gepliver = apply(gs_gepliver, 1, function(i) str_split(i['Top10_Feature_Genes'], ';')[1]) 
gs_list_gepliver = lapply(gs_list_gepliver, function(i) i[[1]]) %>% setNames(gs_gepliver$Cell_Type)
gene_id = data.frame(symbol = c('ACLY', gs_list_gepliver[['B cell']]),
                     entrez = id[match(c('ACLY', gs_list_gepliver[['B cell']]), id$Symbol), ]$ENTREZID)
cor_cell = t(GSE164760.norm[which(rownames(GSE164760.norm) %in% gene_id$symbol), ]) %>% cor()

# Extract ACLY expression from MASH-HCC samples================================
mash_hcc_sample_index = 118:170
acly_id = "47"
acly_exp = GSE164760.exp[acly_id, mash_hcc_sample_index] %>% 
  t() %>% 
  as.data.frame() %>%
  setNames('ACLY')

# Deconvolution=================================================================                         
## Proceed to http://timer.comp-genomics.org/ for immune cell deconvolution using GSE164760_norm.txt.gz with settings species = Human, cancer type = LIHC============
estimation_matrix = read_csv("estimation_matrix.csv")
colnames(estimation_matrix)[-1] = colnames(GSE164760.norm)
estimation_matrix = estimation_matrix[, c(1, match(rownames(acly_exp), colnames(estimation_matrix)))]

# Extract B cell inference
epic_b = subset(estimation_matrix, cell_type == 'B cell_EPIC')[-1] %>% as.numeric()
xcell_b = subset(estimation_matrix, cell_type == 'B cell plasma_XCELL')[-1] %>% as.numeric()
quantiseq_b = subset(estimation_matrix, cell_type == 'B cell_QUANTISEQ')[-1] %>% as.numeric()

##  Alternatively, use immunedeconv package====================================
# EPIC
epic = deconvolute(GSE164760.norm, 'epic')
epic = epic[, c(1, match(rownames(acly_exp), colnames(epic)))]

# XCELL
xcell = deconvolute(GSE164760.norm, 'xcell')
xcell = xcell[, c(1, match(rownames(acly_exp), colnames(xcell)))]

# QUANTISEQ
quantiseq = deconvolute(GSE164760.norm, 'quantiseq')
quantiseq = quantiseq[, c(1, match(rownames(acly_exp), colnames(quantiseq)))]

# Extract B cell inference
epic_b = subset(epic, cell_type == 'B cell')[-1] %>% as.numeric()
xcell_b = subset(xcell, cell_type == 'B cell plasma')[-1] %>% as.numeric()
quantiseq_b = subset(quantiseq, cell_type == 'B cell')[-1] %>% as.numeric()

# Correlate B cell estimate with ACLY expression===============================================
cor.test(epic_b, scale(acly_exp$ACLY))
cor.test(xcell_b, scale(acly_exp$ACLY))
cor.test(quantiseq_b, scale(acly_exp$ACLY))
