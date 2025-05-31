# This script performs differential gene expression analysis using mouse-derived RNA-seq featureCounts data and outputs time-modified result from interaction test and timepoint specific results.
# Output from this script corresponds to Figure 4A, Extended Figure 7C, D, E 
library(dplyr)
library(DESeq2)

# Get file path for counts data
file_path = list.files("/path/to/featureCounts") # modify to folder containing the featureCounts data

# Read
dat = lapply(file_path, function(x) {
  
  fc = read.delim(x, header = T)
  rownames(fc) = fc[, 1]
  fc = fc[, 2, drop = FALSE]
  return(fc)
  
}) %>% do.call('cbind', .)

# Specify conditions
coldata_fc = data.frame(condition = c(rep('WT', 9), rep('KO', 12)), 
                        time = c(rep("Early",5), rep("Late",4), rep("Early",7), rep("Late",5)) 
                       ) # modify according to the sample ordering of your counts matrix

# Overall differential gene expression adjusted for time
deseq2_obj = DESeqDataSetFromMatrix(countData = dat,
                                    colData = coldata_fc,
                                    design = ~ time + condition)
deseq2_obj$condition = relevel(deseq2_obj$condition, ref = "WT")
deseq2_obj = DESeq(deseq2_obj)
res = results(deseq2_obj, test = 'Wald')

# Timepoint specific result and time interaction model
deseq2_obj = DESeqDataSetFromMatrix(countData = dat,
                                    colData = coldata_fc,
                                    design = ~ condition + time + condition:time)
deseq2_obj$condition = relevel(deseq2_obj$condition, ref = "WT")
deseq2_obj$time = relevel(deseq2_obj$time, ref = "Early")
deseq2_obj = DESeq(deseq2_obj)

# Interaction result
res_interaction = results(deseq2_obj, name = "conditionKO.timeLate")

# Early time point result
res_early = results(deseq2_obj, 
                    name = "condition_KO_vs_WT", test = 'Wald'
                    )

# Late time point result
res_late = results(deseq2_obj, 
                   contrast = list(
                     c("condition_KO_vs_WT", 
                       "conditionKO.timeLate")
                     ), test = 'Wald'
                   )
