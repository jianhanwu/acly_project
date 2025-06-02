# This script performs differential gene expression associated with ACLY expression in MASH-HCC using CEL files obtained from GSE164760
# Output from this script corresponds to Extended Figure 12E, F, G
library(dplyr)
library(affy)
library(magrittr)
library(tidyverse)
library(limma)
library(metapredict)
library(hgu219hsentrezgcdf)
library(splines)
library(statmod)

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

# Extract ACLY expression from MASH-HCC samples
mash_hcc_sample_index = 118:170
acly_id = "47"
acly_exp = GSE164760.exp[acly_id, mash_hcc_sample_index] %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  setNames(c('Sample', 'ACLY')) %>%
  arrange(ACLY) %>% 
  mutate(Tertile = c(rep("Low",18), rep('Mid',17), rep('High',18))) %>%
  arrange(Sample) %>% 
  column_to_rownames("Sample")

# Linear model==================================================
# Design matrix
design_acly_linear = model.matrix(~acly_exp$ACLY)

# Fit model
fit_acly_linear <- lmFit(GSE164760.exp[, mash_hcc_sample_index], design_acly_linear)
fit_acly_linear = eBayes(fit_acly_linear)

# Extract genes significantly associated with ACLY expression
result_acly_linear = topTable(fit_acly_linear, coef = 2, n=Inf, sort.by="none") %>%
  rownames_to_column("ENTREZID") %>%
  subset(adj.P.Val <= 0.05)

# PCA
acly_deg_pc = GSE164760.exp[result_acly_linear$ENTREZID, mash_hcc_sample_index] %>%
  t() %>% 
  prcomp(scale = TRUE, center = TRUE)
summary(acly_deg_pc) 

# Nonlinear model====================================================
# Design matrix
acly_spline = ns(acly_exp$ACLY, df= 5)
design_acly = model.matrix(~acly_spline)

# Fit model
fit_acly_spline = lmFit(GSE164760.exp[, mash_hcc_sample_index], design = design_acly)
fit_acly_spline = eBayes(fit_acly_spline, robust = TRUE)

# Extract genes significantly associated with ACLY expression
result_acly_spline = topTable(fit_acly_spline, coef = 2:ncol(design_acly), n = Inf, sort.by = "none") %>%
  rownames_to_column("ENTREZID") %>%
  subset(adj.P.Val <= 0.05)

# PCA
acly_deg_pc = GSE164760.exp[result_acly_spline$ENTREZID, mash_hcc_sample_index] %>%
  t() %>% 
  prcomp(scale = TRUE, center = TRUE)
summary(acly_deg_pc) # nonlinear model more efficiently captures variance (initial PCs explain more variance than equivalent PCs from linear model)
acly_deg_pc = acly_deg_pc$x %>% 
  cbind(Tertile = acly_exp$Tertile) %>% 
  as.data.frame() %>% 
  type.convert(as.is = TRUE)

# Factor model=============================================
# Design matrix (Low vs. High)
acly_bin = factor(c(rep(0,18), rep(1,18)), labels=c("Low","High"))
design_acly_bin = model.matrix(~acly_bin)

# Get sample index matching ACLY bin
bin_sample = acly_exp %>%
  arrange(ACLY, decreasing = F) %>%
  { bind_rows(slice_head(., n = 18), slice_tail(., n = 18)) }
bin_index = match(rownames(bin_sample), colnames(GSE164760.exp))

# Fit model
lmfit_aclybin <- lmFit(GSE164760.exp[, bin_index], design_acly_bin)
contrast_aclybin = makeContrasts(- acly_binHigh,
                                 levels = design_acly_bin)
lmfit.cont_aclybin <- contrasts.fit(lmfit_aclybin, contrast_aclybin)
lmfit_aclybin.ebayes <- eBayes(lmfit_aclybin)

# Extract genes differentially expressed in ACLY high vs. low groups
res_lmfit_aclybin.cont = decideTests(lmfit_aclybin.ebayes, p.value = 0.05)
res_lmfit_aclybin = topTable(lmfit_aclybin.ebayes, coef=2, n=Inf) %>%
  rownames_to_column("ENTREZID")

# CXCL13 expression
cxcl13_exp = GSE164760.exp['10563', bin_index] %>%
  t() %>% as.data.frame() %>% 
  setNames('CXCL13') %>%
  cbind(Bin = bin_sample$Tertile)
lmfit_aclybin.cont_p = lmfit_aclybin.ebayes$p.value %>%
  as.data.frame() %>%
  rownames_to_column("ENTREZID") %>%
  filter(ENTREZID == '10563') %>%
  select(contains('acly_bin'))
