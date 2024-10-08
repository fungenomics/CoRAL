# load libraries and arguments 
library(scPred)
library(Seurat)
library(tidyverse)
library(WGCNA)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
query_path = args[1]
model_path = args[2]
pred_path = args[3]
threads = as.numeric(args[4])
model_type = args[5]
threshold = as.numeric(args[6])
# get path for other output
out_path = dirname(pred_path)

#--------------- Data -------------------

# read query matrix 
message('@ READ QUERY')
query = data.table::fread(query_path, nThread=threads, header=T, data.table=F) %>%
        column_to_rownames('V1') 
message('@ DONE')

# load model 
message('@ LOAD MODEL')
load(model_path)
message('@ DONE')

# transpose and create seurat object 
query = transposeBigData(query)
query = CreateSeuratObject(query)

# normalize query 
query = query %>% 
  NormalizeData() 

#----------- Predict scPred -------------

# predict cells 
message('@ PREDICT LABELS')
query = scPredict(query,
                  scpred,
                  threshold = threshold)
message('@ DONE')

head(colnames(query))
head(query$scpred_prediction)

# scPred chnages - to _minus --> chnage back before saving 
query$scpred_prediction = gsub("_minus", "-", query$scpred_prediction)
pred_labs = data.frame(cell = colnames(query),
                       scPred = query$scpred_prediction)
colnames(pred_labs)[2] = paste0('scPred_',model_type)
# write prediction 
data.table::fwrite(pred_labs, file = pred_path)

# save probbability matrix 
prob_mat = query@meta.data %>% select(starts_with('scpred'))
colnames(prob_mat) = str_remove(colnames(prob_mat), 'scpred_')
colnames(prob_mat) = gsub("_minus", "-", colnames(prob_mat))
# Remove the column that are not prob informative 
prob_mat = prob_mat[,-which(colnames(prob_mat) %in% c("max","prediction","no_rejection"))]
## Normalization to use in CAWPE
prob_mat[prob_mat < 0] <- 0
prob_mat <- apply(prob_mat,1,function(x){
  x / sum(x)
}) %>% t()
prob_mat = prob_mat %>% as.data.frame() %>% rownames_to_column('cell')
colnames(prob_mat)[1] = ""

# write probability matrix 
data.table::fwrite(prob_mat, file = paste0(out_path, '/scPred_',model_type,'_pred_score.csv'))

#----------------------------------------
