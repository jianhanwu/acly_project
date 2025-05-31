# This script prepares TPM normalized expression matrix for use with TIMER2.0 immune deconvolution
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

# NOTE: Deconvolution analyses were performed using TIMER2.0 hosted at http://timer.cistrome.org/ at the time of analysis. However, the server seems to be dysfunct at the time of publication. 
# Nevertheless, results can be reproduced by using the immunedeconv package as shown below. 
# Clean matrix for deconvolution
exp_dat = aclyko_tpm[, -which(colnames(aclyko_tpm) %in% c('ENTREZID', 'gene_symbol'))]
rownames(exp_dat) = aclyko_tpm[, 'gene_symbol']

# mMCP-counter
mmcp = deconvolute_mouse(exp_dat, "mmcp_counter")

# Correlation with ACLY expression
fc_tpm_mcp = data.frame(ACLY = as.numeric(subset(aclyko_tpm, gene_symbol == 'Acly')[-c(1,2)]),  
                        B_cell = as.numeric(mmcp[4, -1]), 
                        condition = c(rep('KO', 12), rep('WT', 9)) #modify based on sample ordering
                        )
cor.test(fc_tpm_mcp$ACLY, fc_tpm_mcp$B_cell)
