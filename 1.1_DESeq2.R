# This script performs differential gene expression analysis using mouse-derived RNA-seq featureCounts data and outputs time-modified result from interaction test and timepoint specific results.
library(dplyr)
library(DESeq2)

# Get file path for counts data
file_path = list.files("/path/to/featureCounts") # modify to folder containing the featureCounts data

# Read
dat = sapply(file_path, function(x) {
  
  fc = read.delim(x, header = T)
  rownames(fc) = fc[, 1]
  fc = fc[, -1]
  
}) %>% do.call('cbind', .)

# Specify conditions
coldata_fc = data.frame(condition = c(rep('WT', 9), rep('KO', 12)), 
                        time = c(rep("Early",5), rep("Late",4), rep("Early",7), rep("Late",5)) 
                       ) # modify according to the sample ordering of your counts matrix

# DESeq2
deseq2_obj = DESeqDataSetFromMatrix(countData = dat,
                                     colData = coldata_fc,
                                     design = ~ condition*time)

# Set reference level
deseq2_obj$condition = relevel(deseq2_obj$condition, ref = "WT")
deseq2_obj$Time = relevel(deseq2_obj$Time, ref = "Early")

# Run DESeq2
deseq2_obj = DESeq(deseq2_obj, test="LRT", reduced = ~ condition + time)

# Interaction result
res_interaction = results(deseq2_obj, name = "conditionKO.timeLate")

# Early time point result
res_early = results(deseq2_obj, name="condition_KO_vs_WT")
res_early = lfcShrink(deseq2_obj, contrast="condition_KO_vs_WT", res = res_early, type = "normal")

# Late time point result
res_late = results(deseq2_obj, 
                   contrast = list(
                     c("condition_KO_vs_WT", 
                     "conditionKO.timeLate")
                   )
                  )
res_late = lfcShrink(deseq2_obj, 
                     contrast = list(
                       c("condition_KO_vs_WT", 
                         "conditionKO.timeLate")
                     )
                     res = res_late, type = "normal")
