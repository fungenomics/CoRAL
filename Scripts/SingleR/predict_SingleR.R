# load libraries
library(tidyverse)
library(SingleR)
library(SingleCellExperiment)
library(WGCNA)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
sample_path = args[1]
model_path = args[2]
pred_path = args[3]
threads = as.numeric(args[4])

# get path for other output
out_path = dirname(pred_path)

#--------------- Data -------------------

# read query matrix and transpose 
message('@ READ QUERY')
query = data.table::fread(sample_path, nThread=threads, header=T, data.table=F) %>%
        column_to_rownames('V1') 
message('@ DONE')

# load model 
message('@ LOAD MODEL')
load(model_path)
message('@ DONE')

# Make SingleCellExperiment object from query (transpose query first)
query = transposeBigData(query)
query = SingleCellExperiment(assays = list(counts = query))

# Log normalize query counts 
message('@ NORMALIZE QUERY')
query = scuttle::logNormCounts(query)
message('@ DONE')

#----------- Predict SingleR ------------

# specify parallelization configuration depending on number of threads
if(threads > 1){
  
  bpparam <- BiocParallel::MulticoreParam(workers = threads)
  # parallel <- TRUE
  
} else {
  
  bpparam <- BiocParallel::SerialParam()
  # parallel <- FALSE
  
}

# predict labels with the fine tune method
message('@ PREDICT LABELS')
pred = classifySingleR(query,
                       singler,
                       assay.type = "logcounts",
                       BPPARAM = bpparam)
message('@ DONE')

# head(pred)

pred_labs = data.frame(cell = rownames(pred),
	                     SingleR = pred$pruned.labels)
#Pruned labels are reported in the pruned.labels field where
#low-quality assignments are replaced with NA.
pred_labs$SingleR[is.na(pred_labs$SingleR)] <- 'Unknown'
# write prediction 
data.table::fwrite(pred_labs, file = pred_path)

# save pred output
## Save the original correlation matrix
save(pred, file = paste0(out_path, '/SingleR_output.Rdata'))

#----------------------------------------
# The scores matrix has several caveats associated with its interpretation.
# Only the pre-tuned scores are stored in this matrix,
#as scores after fine-tuning are not comparable across all labels.
# This means that the label with the highest score for a cell may not be the
#cell’s final label if fine-tuning is applied. (https://bioconductor.org/books/devel/SingleRBook/annotation-diagnostics.html)
# Because of this we use the pre-fine-tunning method to compute the probabilities for CAWPE 
## Normalize the correlation matrix
prob_matrix = pred$scores %>% `rownames<-`(rownames(pred))
prob_matrix[prob_matrix < 0] <- 0
prob_matrix <- apply(prob_matrix,1,function(x){
  x / sum(x)
}) %>% t()

prob_matrix = prob_matrix %>% as.data.frame() %>% rownames_to_column('cell')
colnames(prob_matrix)[1] = ""

data.table::fwrite(prob_matrix,
                   file = paste0(out_path, '/SingleR_pred_score.csv'))