# This script evaluates ACLY expression in normal and malignant hepatocytes using snRNA-seq data from GSE189175
# Output from this script corresponds to Figure 5O
library(data.table)
library(dplyr)
library(Seurat)
library(glmGamPoi)

#1. Load expression matrix and create Seurat objec===========================
expression_matrix = Read10X(data.dir = './GSE189175/', cell.column = 1)
seurat_object <- CreateSeuratObject(counts = expression_matrix, project = 'GSE189175')
rm_genes = which(is.na(features$V2))
seurat_object = subset(seurat_object, features = features[-rm_genes,]$V2)

#2. QC=============
seurat_object@active.assay = 'RNA'
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-") #0 mitochondrial DNA as it is nuclei seq
seurat_object[["pct.alb"]] <- PercentageFeatureSet(seurat_object, features = "ALB") 

#3. SCTransform====================
vars.to.regress = NULL
seurat_object = SCTransform(seurat_object, conserve.memory = TRUE, vars.to.regress = vars.to.regress, verbose = F)

#4. Dimension Reduction and Cluster==============================================
seurat_object@active.assay = 'SCT'
seurat_object = RunPCA(seurat_object, npcs = 30, verbose = FALSE) 
seurat_object <- RunUMAP(seurat_object, dims = 1:30, verbose = FALSE)

# Find clusters
seurat_object <- FindNeighbors(seurat_object, dims = 1:30, verbose = TRUE)
seurat_object <- FindClusters(seurat_object, resolution = c(0.2, 0.5, 1, 1.5, 2), verbose = TRUE)

#5. Assign cell types using sc-type================================================
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# prepare gene sets
gs_gepliver = read.csv('./gepliver_cell_annot.csv')
gs_list_gepliver = apply(gs_gepliver, 1, function(i) str_split(i['Top10_Feature_Genes'], ';')[1]) 
gs_list_gepliver = lapply(gs_list_gepliver, function(i) i[[1]]) %>% setNames(gs_gepliver$Cell_Type)
gs_list_gepliver = list(gs_positive = gs_list_gepliver)

# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 = isFALSE('counts' %in% names(attributes(seurat_object[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
scRNAseqData_scaled = if (seurat_package_v5) as.matrix(seurat_object[["SCT"]]$scale.data) else as.matrix(seurat_object[["SCT"]]@scale.data)

# run ScType
es.max_gepliver = sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list_gepliver$gs_positive)
es.max_gepliver_major = es.max_gepliver[2:17, ]
  
# merge by cluster
es.max = es.max_gepliver_major
cL_resutls_gep <- do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_object@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores_gep <- cL_resutls_gep %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores_gep$type[as.numeric(as.character(sctype_scores_gep$scores)) < sctype_scores_gep$ncells/4] <- "Unknown"
print(sctype_scores_gep[,1:3])

# Assign cell type to cells in seurat object
seurat_object@meta.data$gep_classification = ""
for(j in unique(sctype_scores_gep$cluster)){
  cl_type = sctype_scores_gep[sctype_scores_gep$cluster==j,]; 
  seurat_object@meta.data$gep_classification[seurat_object@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

#6. Evaluate ACLY expression in cell types===============================================
exp = GetAssayData(seurat_object, slot = 'data')

gene = 'ACLY'
acly_exp = exp[which(rownames(exp) == gene), ]
all(names(acly_exp) == rownames(seurat_object@meta.data)) # TRUE
acly_exp = cbind(seurat_object@meta.data[, !grepl('SCT_snn', colnames(seurat_object@meta.data))], gene = acly_exp)
acly_exp$tissue = ifelse(grepl('hcc', rownames(acly_exp)), 'MASH-HCC', 'MASH-HCC Adj')

# Frequency of gene expression per cell types
nashhcc_freq = acly_exp %>% group_by(tissue, gep_classification) %>%
  dplyr::summarise(n_cell = length(gene), freq_exp = length(which(gene != 0))/ length(gene))

# Proportion test, MASH-HCC vs MASH-HCC adjacent tissue
dat = nashhcc_freq
colnames(dat)[2] = 'Cell_Type'
dat$n_cell_exp = dat$n_cell*dat$freq_exp
dat$n_cell_not_exp = dat$n_cell - dat$n_cell*dat$freq_exp

prop_test = lapply(names(which(table(dat$Cell_Type) == 2)), function(i) {
  
  message(i)
  if (i != 'Hepa_Malignant') {
    # Compare proportion of cells expressing ACLY between HCC and HCC adjacent tissue
    p_res = prop.test(x = subset(dat, Cell_Type == i)[1:2,]$n_cell_exp, 
                      n = subset(dat, Cell_Type == i)[1:2,]$n_cell)
    prop_dat = c(p_res$statistic, p = p_res$p.value, Cell_Type = i)
  
  } else {
    
    # Compare proportion of cells expressing ACLY between normal and malignant hepatocytes
    dat_hepa = subset(dat, grepl('Hepa', Cell_Type))
    p_res = prop.test(x = dat_hepa[1:2,]$n_cell_exp, 
                      n = dat_hepa[1:2,]$n_cell)
    prop_dat = c(p_res$statistic, p = p_res$p.value, Cell_Type = i)
  }
  
}) %>% do.call('rbind',.)
prop_test = as.data.frame(prop_test)
prop_test = type.convert(prop_test)
prop_test$Condition = ifelse(prop_test$Cell_Type == 'Hepa_Malignant', 'Hepa_Malignant vs. Hepa_Normal', 'HCC vs. Adj')

annot = data.frame(Cell_Type = prop_test$Cell_Type, 
                   freq_exp = subset(dat, tissue == 'MASH-HCC')[match(prop_test$Cell_Type, subset(dat, tissue == 'MASH-HCC')$Cell_Type),]$freq_exp,
                   lab = paste0('p = ', signif(prop_test$p, 1)), 
                   tissue = 'MASH-HCC', 
                   Condition = prop_test$Condition 
                   )
                          
# Extract hepatocyte specific result
dat_hep = subset(dat, grepl('Hepa', Cell_Type))
annot_hep = subset(annot, grepl('Hepa', Cell_Type))
