library(tidyverse)
library(Seurat)
library(singleCellNet)
library(tidyverse)
library(data.table)
library(WGCNA)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)
sample_path = args[1]
model_path = args[2]
pred_path = args[3]
threads = as.numeric(args[4])
threshold = as.numeric(args[5])
# get path for other output
out_path = dirname(pred_path)

#--------------- Data -------------------

# read query matrix and transpose 
message('@ READ QUERY')
query = data.table::fread(sample_path, nThread=threads, header=T, data.table=F) %>%
        column_to_rownames('V1') 
message('@ DONE')

# Get cell names
#cellnames = row.names(query)

# load model 
message('@ LOAD MODEL')
load(model_path)
message('@ DONE')

# Transpose query 
query = transposeBigData(query, blocksize = 10000)

#----------- Predict singleCellNet ------------

# predict labels 
message('@ PREDICT LABELS')

# No generating random cells
pred = scn_predict(class_info[['cnProc']], query, nrand = 0)
message('@ DONE')

# classify cells , I do the same as what they do in get_cate, 
# assign the cell to the cell-type with the highest score
# if the score is lower than certain threshold, assign it as Unknown (rand)
df <- t(pred) %>% 
  apply(1,function(x){
    data.frame(label = names(x)[which.max(x)],
               score = max(x))
    }) %>% 
  do.call(what = "rbind",.) %>% 
  mutate(label = ifelse(label == "rand","Unknown",label)) #Convert rand to Unknown
#It is assisgned rand when the category with the max score is rand, then the rand is convert to Unknown
# finally for those assignation when the max value is lower than a certain threshold (0.5 by default)
# the label is also converted to Unknown
df$label[df$score < threshold] <- "Unknown" 

pred_labs = data.frame(cell = rownames(df),
	                     singleCellNet = df$label)
rm(df)
# write prediction 
data.table::fwrite(pred_labs,
                   file = pred_path)

# The prob matrix includes the rand category, since we are planning to
# use it in CAWPE, we should remove the rand category and re-normalized.
# removing rand column
prob_matrix = t(pred) %>% .[,-which(colnames(.) == "rand")]
prob_matrix[prob_matrix < 0] <- 0
prob_matrix <- apply(prob_matrix,1,function(x){
  x / sum(x)
}) %>% t()

prob_matrix = prob_matrix %>% as.data.frame() %>% rownames_to_column('cell')
colnames(prob_matrix)[1] = ""


data.table::fwrite(prob_matrix, file = paste0(out_path, '/singleCellNet_pred_score.csv'))

#----------------------------------------
