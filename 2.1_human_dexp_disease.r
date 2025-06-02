# This script performs differential gene expression analysis between disease states using CEL files obtained from GSE164760
# Output from this script corresponds to Figure 5N, Extended Figure 12A, B, C
library(dplyr)
library(affy)
library(magrittr)
library(tidyverse)
library(limma)
BiocManager::install('metapredict', site_repository = 'https://hugheylab.github.io/drat/')
library(metapredict)
installCustomCdfPackages("hgu219hsentrezgcdf", ver = 25)
library(hgu219hsentrezgcdf)
library(ggVennDiagram)

# Load gene ID
id = read.delim('hgu219_entrez_symbol.txt', header = T)

# Get file path for RMA data
dir_path = '/path/to/GSE164760/'

# Read
GSE164760 = ReadAffy(celfile.path = dir_path, cdfname="hgu219hsentrezg")

# Create epression file
norm.GSE164760 = rma(GSE164760)
GSE164760.exp = as.data.frame(exprs(norm.GSE164760))
rownames(GSE164760.exp) = gsub('_at', '', rownames(GSE164760.exp))

# Design matrix
tissue = factor(c(rep(1,6), rep(2,8), rep(3,74), rep(4,29), rep(5,53)),
                labels=c("Healthy","Cirrhosis","MASH","NTA", "HCC"))
design_GSE164760 = model.matrix(~0 + tissue)
colnames(design_GSE164760) = c("Healthy","Cirrhosis","MASH","NTA", "HCC")

# Fit model
lmfit <- lmFit(GSE164760.exp, design_GSE164760)
contrast_GSE164760 = makeContrasts(HCC-MASH, 
                                   HCC-NTA, 
                                   NTA-MASH, 
                                   HCC-Healthy, 
                                   HCC-Cirrhosis, 
                                   MASH-Healthy, 
                                   NTA-Healthy, 
                                   Cirrhosis-Healthy,
                                   levels = design_GSE164760)
lmfit.cont <- contrasts.fit(lmfit, contrast_GSE164760)
lmfit.cont.ebayes <- eBayes(lmfit.cont)

# Pairwise contrasts
# HCC vs. MASH
res_lmfit.cont_hccmash = topTable(lmfit.cont.ebayes, coef=1, n=Inf) %>% 
  rownames_to_column("ENTREZID")
# HCC vs. NTA
res_lmfit.cont_hccnta = topTable(lmfit.cont.ebayes, coef=2, n=Inf) %>% 
  rownames_to_column("ENTREZID")
# HCC vs. Healthy
res_lmfit.cont_hcchealthy = topTable(lmfit.cont.ebayes, coef=4, n=Inf) %>% 
  rownames_to_column("ENTREZID") 
# HCC vs. Cirrhosis
res_lmfit.cont_hcccirrhosis = topTable(lmfit.cont.ebayes, coef=5, n=Inf) %>% 
  rownames_to_column("ENTREZID")

# Find genes differentially expressed across all pairwise MASH-HCC comparisons
get_sig_genes <- function(data, direction = "UP") {
  if (direction == "UP") {
    subset(data, adj.P.Val <= 0.05 & logFC > 0)[, 1]
  } else {
    subset(data, adj.P.Val <= 0.05 & logFC < 0)[, 1]
  }
}
list_limma <- list(
  'HCC-NASH_UP'       = get_sig_genes(res_lmfit.cont_hccmash, "UP"),
  'HCC-NTA_UP'        = get_sig_genes(res_lmfit.cont_hccnta, "UP"),
  'HCC-Healthy_UP'    = get_sig_genes(res_lmfit.cont_hcchealthy, "UP"),
  'HCC-Cirrhosis_UP'  = get_sig_genes(res_lmfit.cont_hcccirrhosis, "UP"),
  'HCC-NASH_DN'       = get_sig_genes(res_lmfit.cont_hccmash, "DN"),
  'HCC-NTA_DN'        = get_sig_genes(res_lmfit.cont_hccnta, "DN"),
  'HCC-Healthy_DN'    = get_sig_genes(res_lmfit.cont_hcchealthy, "DN"),
  'HCC-Cirrhosis_DN'  = get_sig_genes(res_lmfit.cont_hcccirrhosis, "DN")
)
common_up = process_region_data(Venn(list_limma[c(1:4)]))
common_up_genes = common_up$item[[15]]
common_dn = process_region_data(Venn(list_limma[-c(1:4)]))
common_dn_genes = common_dn$item[[15]]

# Extract ACLY, ACACA, ACACB, ACSS2, FASN differential expression result
entrezid_symbol = c(ACLY = 47, ACACA = 31, ACACB = 32, ACSS2 = 55902, FASN = 2194)
intersect(common_up_genes, entrezid_symbol)
intersect(common_dn_genes, entrezid_symbol)
select_dexp_res = lmfit.cont.ebayes$p.value %>%
  as.data.frame() %>%
  rownames_to_column("ENTREZID") %>%
  filter(ENTREZID %in% as.numeric(entrezid_symbol)) 


