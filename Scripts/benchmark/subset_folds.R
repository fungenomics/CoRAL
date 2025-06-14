# load libraries and arguments 
library(rBayesianOptimization)
library(tidyverse)

initial.options = commandArgs(trailingOnly = FALSE)
file.arg.name = "--file="
script.name = sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)]) 

source(paste0(dirname(dirname(script.name)), "/Functions/functions.R"))

args = commandArgs(trailingOnly = TRUE)

print(paste0("@ Using seed ",as.character(args[10])))
set.seed(as.numeric(args[10]))

ref_path = args[1]
lab_path = args[2]
out_path = args[3]

threads = as.numeric(args[4])
if(is.na(threads)){
  stop("The number threads specified is not a numeric value")
}

n_folds = as.numeric(args[5])
if(is.na(n_folds)){
  stop("The number of folds specified is not a numeric value")
}

min_cells = as.numeric(args[6])
if(is.na(min_cells)){
  stop("The minimum number of cells specified is not a numeric value")
}

downsample_value = as.numeric(args[7]) 
if(is.na(downsample_value)){
  stop("The downsample value specified is not a numeric value")
}
print(class(downsample_value))
downsample_stratified = as.logical(args[8])
if(is.na(downsample_stratified)){
  stop("The downsample stratified specified is not a logical value")
}

#downsample_stratified = if(downsample_stratified) "label" else NULL

# ontology_path = args[9]
# ontology_columns = strsplit(args[10], split = ' ')[[1]]

batch_path = args[9]
if(batch_path == 'None'){
  batch_path = NULL
}

print(batch_path)
print(downsample_stratified)
print(class(downsample_stratified))
# print(ontology_path)
# print(ontology_columns)

#--------------- Data -------------------
# read reference matrix 
message('@ READ REF')
tmp <- get_data_reference(ref_path = ref_path,
                          lab_path = lab_path,
                          batch_path = batch_path)
ref     <- tmp$exp
labels  <- tmp$lab
rm(tmp)
message('@ DONE')

# downsample 
if(downsample_value != 0){
  labels = downsample(labels, downsample_stratified, downsample_value)
}

# remove small clusters 
if(min_cells > 0){
  labels = remove_small_clusters(labels, min_cells)
}

ref = ref[rownames(labels),]

# save downsampled lables 
save.df <- data.frame(cells= rownames(labels), 
                      labels)

colnames(save.df)[1] <- ""

data.table::fwrite(save.df,
                   file = paste0(out_path,'/downsampled_reference_labels.csv'),
                   col.names = T,
                   row.names=F,
                   sep = ",")
rm(save.df)

# check if cell names are in the same order in labels and ref
order = all(as.character(rownames(labels)) == as.character(rownames(ref)))

# throw error if order is not the same 
if(!order){
    stop("@ Order of cells in reference and labels do not match")
}

# ref[1:10, 1:10]
# head(labels)

# create n folds 
folds = KFold(labels$label, 
              nfolds = n_folds, 
              stratified = T, 
              seed = 1234)
head(folds)

# Loop through folds and save training and testing data sets 
for (i in 1:n_folds){
  message(paste0('@ SAVING FOLD ', i))

  # subset test fold
  message('subset test fold')
  test = ref[folds[[i]], ,drop=F]
  test = test %>% rownames_to_column("cell")
  colnames(test)[1] = ""

  # subset true test labels 
  message('subset true test labels')
  test_labels = labels[folds[[i]], ,drop=F]
  test_labels = test_labels %>% rownames_to_column("cell")
  colnames(test_labels)[1] = ""
  
  # subset training data 
  message('subset true test labels')
  train = ref[-folds[[i]], ,drop=F]
  train = train %>% rownames_to_column("cell")
  colnames(train)[1] = "" 
   
  # subset labels for training data
  message('@ subset labels for training data')
  train_labels = labels[-folds[[i]], ,drop=F]
  train_labels = train_labels %>% rownames_to_column("cell")
  colnames(train_labels)[1] = ""
  
  # check if you have enough cells per label in each fold
  min_cell_per_fold_ct <- 10
  if(!all(table(train_labels$label) >= min_cell_per_fold_ct)){
    stop(paste0("In fold ",i, " not all the training labels have more than ", min_cell_per_fold_ct, " cells"))
  }
  # save csv files 
  data.table::fwrite(test, paste0(out_path, '/fold', i, '/test.csv'))
  data.table::fwrite(test_labels, paste0(out_path, '/fold', i, '/test_labels.csv'))
  data.table::fwrite(train, paste0(out_path, '/fold', i, '/train.csv'))
  data.table::fwrite(train_labels, paste0(out_path, '/fold', i, '/train_labels.csv'))
}
