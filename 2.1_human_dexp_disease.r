# This script performs differential gene expression analysis between disease states using CEL files obtained from GSE164760
# Output from this script corresponds to Figure 5N, Extended Figure 12A, B, C
library(dplyr)
library(affy)
library(magrittr)
library(tidyverse)
library(limma)
library(GEOquery)
BiocManager::install('metapredict', site_repository = 'https://hugheylab.github.io/drat/')
library(metapredict)
installCustomCdfPackages("hgu219hsentrezgcdf", ver = 25)
library(hgu219hsentrezgcdf)

# Load probe ID
id = read.delim('GSE164760_probe_id.txt', header = F)
colnames(id) = c('Probe', 'ENTREZID', 'Symbol')

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
