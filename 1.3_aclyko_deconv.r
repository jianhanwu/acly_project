# This script prepares TPM normalized expression matrix for immune cell deconvolution
# Output from this script corresponds to Figure 4D
library(dplyr)
library(org.Mm.eg.db)
library(mMCPcounter)
library(magrittr)
library(tidyverse)

# Get file path for counts data
file_path = list.files("/path/to/featureCounts") # modify to folder containing the featureCounts data

# Read
dat = lapply(file_path, function(x) {
  
  fc = read.delim(x, header = T)
  rownames(fc) = fc[, 1]
  fc = fc[, 2, drop = FALSE]
  return(fc)
  
}) %>% do.call('cbind', .)

# Read gene length data 
fc_length = readRDS('path/to/feature_lengths/fc_length.RDS')

# Convert to TPM
aclyko_tpm = do.call("cbind", lapply(1:21, function(x) {
  ratio = dat[,x]/fc_length[,2]
  ratioSum = sum(ratio)
  tpm = 1e6*ratio/ratioSum
})) %>%
  set_colnames(colnames(dat)) %>%
  set_rownames(rownames(dat)) %>%
  as.data.frame()

# Convert entrez ID to gene symbol
aclyko_tpm = rownames_to_column(aclyko_tpm, 'ENTREZID')
symbol = data.frame(
  "gene_symbol" = mapIds(
    org.Mm.eg.db,
    keys = as.character(fc_length$Geneid),
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
) %>%
  filter(!is.na(gene_symbol)) %>%
  rownames_to_column('ENTREZID')
aclyko_tpm = inner_join(symbol, aclyko_tpm, by = 'ENTREZID')

# Clean matrix for deconvolution input
exp_dat = aclyko_tpm[, -which(colnames(aclyko_tpm) %in% c('ENTREZID', 'gene_symbol'))]
rownames(exp_dat) = aclyko_tpm[, 'gene_symbol']
write.table(exp_dat, './aclyko_tpm.txt', col.names = T, row.names = T, quote = F, sep = '\t')

# Deconvolution=============================================================================
# Proceed to http://timer.comp-genomics.org/ for immune cell deconvolution with settings species = Mouse, cancer type = LIHC===========
estimation_matrix = read_csv("aclyko_estimation_matrix.csv")

# Extract B cell inference
mmcp_b = subset(estimation_matrix, cell_type == 'B cell_MMCPCOUNTER')[-1] %>% as.numeric()
timer_b = subset(estimation_matrix, cell_type == 'B cell_TIMER')[-1] %>% as.numeric()
cibersort_b = subset(estimation_matrix, cell_type == 'B cell plasma_CIBERSORT-ABS')[-1] %>% as.numeric()
quantiseq_b = subset(estimation_matrix, cell_type == 'B cell_QUANTISEQ')[-1] %>% as.numeric()

# Alternatively use immunedeconv package============================================
# mMCP-counter
mmcp = deconvolute_mouse(exp_dat, "mmcp_counter")
mmcp_b = as.numeric(mmcp[4, -1])

# TIMER
exp_dat_h_ortho <- immunedeconv::mouse_genes_to_human(exp_dat)
timer = deconvolute(exp_dat_h_ortho, 'timer')

# CIBERSORT-ABS
cibersort = deconvolute(exp_dat_h_ortho, 'cibersort_abs')

# QUANTISEQ
quantiseq = deconvolute(exp_dat_h_ortho, 'quantiseq')

# Correlation with ACLY expression====================================================
acly_exp = as.numeric(subset(aclyko_tpm, gene_symbol == 'Acly')[-c(1,2)])
cor.test(mmcp_b, acly_exp)
cor.test(timer_b, acly_exp)
cor.test(cibersort_b, acly_exp)
cor.test(quantiseq_b, acly_exp)
