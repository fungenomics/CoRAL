
library(tidyverse)
library(glue)
initial.options = commandArgs(trailingOnly = FALSE)
file.arg.name = "--file="
script.name = sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)]) 

source(paste0(dirname(script.name), "/Functions/functions.R"))

args = commandArgs(trailingOnly = TRUE)

print(paste0("@ Using seed ",as.character(args[14])))
set.seed(as.numeric(args[14]))

ref_path = args[1]
query_paths = strsplit(args[2], split = ' ')[[1]]
out = args[3]
lab_path = args[5]
reference_name = args[6]
query_names = strsplit(args[7], split = ' ')[[1]]

min_cells = as.numeric(args[8])
if(is.na(min_cells)){
  stop("The minimum number of cells specified is not a numeric value")
}

downsample_value = as.numeric(args[9]) 
if(is.na(downsample_value)){
  stop("The downsample value specified is not a numeric value")
}

downsample_per_class = as.logical(args[10])
if(is.na(downsample_per_class)){
  stop("The downsample stratified specified is not a logical value")
}

names(query_paths) = query_names

batch_path = args[11]
if(batch_path == 'None'){
  batch_path = NULL
}

print(batch_path)

feature_selection_method = strsplit(args[12], split = ' ')[[1]]
if(any(!(feature_selection_method %in% c("intersection","complete")))){
  stop("The method for feature selection should be intersection or complete")
}
print(feature_selection_method)

pipeline_mode = args[13]
print(pipeline_mode)

# threshold overlap genes 
gene_overlap_threshold = as.numeric(args[15])

# Reference gene convertion parameters
convert.genes.ref = as.logical(args[4])
if(convert.genes.ref){
  ref.gene_map <- list(from_species = as.character(args[16]),
                       from_gene = as.character(args[17]),
                       to_species = as.character(args[18]),
                       to_gene = as.character(args[19])
  )
}

if(pipeline_mode == 'annotation'){
  convert.genes.query = as.logical(args[20])
  if(convert.genes.query){
    query.gene_map <- list(from_species = as.character(args[21]),
                           from_gene = as.character(args[22]),
                           to_species = as.character(args[23]),
                           to_gene = as.character(args[24])
    )
  }  
}

# ----- PREPROCESS REFERENCE ----------------------
tmp <- get_data_reference(ref_path = ref_path,
                          lab_path = lab_path,
                          batch_path = batch_path)
data <- list()
data[['ref']] <- tmp$exp
lab           <- tmp$lab
rm(tmp)

# downsample 
if(downsample_value != 0){
  lab = downsample(lab, downsample_per_class, downsample_value)
}

# remove small clusters 
if(min_cells > 0){
  lab = remove_small_clusters(lab, min_cells)
}
# filter reference for donwsampled cells 
data[['ref']] = data[['ref']][rownames(lab),]

# save downsampled lables 
save.df <- data.frame(cells= rownames(lab), 
                      lab)

colnames(save.df)[1] <- ""

data.table::fwrite(save.df,
                   file = paste0(out, '/model/', reference_name, '/downsampled_labels.csv'),
                   col.names = T,
                   row.names=F,
                   sep = ",")
rm(save.df)

# if specified by user, convert reference gene names from mouse to human
if(convert.genes.ref){
  
  message('@ CONVERTING GENE NAMES')

  # include functions and libraries for conversion
  library(Orthology.eg.db)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(WGCNA)

  # convert
  mapped.df = mapfun(genes = colnames(data[['ref']]),
              from = ref.gene_map$from_species,
              from_gene = ref.gene_map$from_gene,
              to = ref.gene_map$to_species,
              to_gene = ref.gene_map$to_gene
              ) %>% dplyr::select(from_genes, to_genes) 
  
  # output list of mouse genes that were not converted
  not_converted = mapped.df %>% filter(is.na(to_genes)) %>% .$from_genes
  
  data.table::fwrite(as.list(not_converted), file = paste0(out, '/model/', reference_name, '/genes_not_converted.csv'), sep = ',')

  # throw error if more than threshold % genes not converted
  threshold = 0.5
  if(length(not_converted) > threshold*length(colnames(data[['ref']]))){
    stop(paste0("@ More than ",threshold*100,"% of mouse genes in reference could not be converted to human"))
  }

  # modify reference matrix to contain converted genes
  data[['ref']] = data[['ref']] %>%
    transposeBigData() %>%
    rownames_to_column('from_genes') %>%
    inner_join(mapped.df %>% filter(!is.na(to_genes)), 
               by = 'from_genes') %>%
    dplyr::select(-from_genes) %>%
    column_to_rownames('to_genes') %>%
    transposeBigData() 
  
}

# ----- QUERY --------------------------------
if(pipeline_mode == "annotation"){
  # read query 
  for(i in 1:length(query_paths)){
  
    print(query_paths[i])
    tmp = get_data_query(query_path = query_paths[i])
    query = names(query_paths)[i]
  
    print(query)
    data[[query]] = tmp
    
    # if specified by user, convert reference gene names from mouse to human
    if(convert.genes.query){
      
      message('@ CONVERTING GENE NAMES')
      
      # include functions and libraries for conversion
      library(Orthology.eg.db)
      library(org.Mm.eg.db)
      library(org.Hs.eg.db)
      library(WGCNA)
      
      # convert
      mapped.df = mapfun(genes = colnames(data[[query]]),
                         from = query.gene_map$from_species,
                         from_gene = query.gene_map$from_gene,
                         to = query.gene_map$to_species,
                         to_gene = query.gene_map$to_gene
      ) %>% dplyr::select(from_genes, to_genes) 
      
      # output list of mouse genes that were not converted
      not_converted = mapped.df %>% filter(is.na(to_genes)) %>% .$from_genes
      
      data.table::fwrite(as.list(not_converted), file = paste0(out, '/', query, '/', reference_name, '/genes_not_converted.csv'), sep = ',')
      # throw error if more than threshold % genes not converted
      threshold = 0.5
      if(length(not_converted) > threshold*length(colnames(data[[query]]))){
        stop(paste0("@ More than ",threshold*100,"% of mouse genes in reference could not be converted to human"))
      }
      
      # modify reference matrix to contain converted genes
      data[[query]] = data[[query]] %>%
        transposeBigData() %>%
        rownames_to_column('from_genes') %>%
        inner_join(mapped.df %>% filter(!is.na(to_genes)), 
                   by = 'from_genes') %>%
        dplyr::select(-from_genes) %>%
        column_to_rownames('to_genes') %>%
        transposeBigData() 
      
    }
  }
}

# ----- GENE SELECTION ------------------------
for(ft_mth in feature_selection_method){
  if(ft_mth == "intersection"){
    ### Takes the intersection between the query and reference feature space
    # get genes for each data frame (colnames)
    genes = lapply(data, function(x){(colnames(x))})
    # reduce set of genes to the intersect 
    common_genes = Reduce(intersect,genes)
    print(paste0('@Found ', length(common_genes), ' in common'))
    
    # throw error if number of common genes below % threshold of genes in any of provided datasets (ref or query) 
    frac = lapply(genes, function(x){length(common_genes)/length(x)})
    if(any(frac < gene_overlap_threshold)){
      names(frac) = names(data)
      print(frac)
      stop(paste0("@ In at least one provided dataset (ref or query), less than ",gene_overlap_threshold*100,"% of genes appear in common gene set. See above for the fraction of genes from each dataset appearing in common gene set (note: samples with few genes will have higher fractions)"))
    }
    # save common genes 
    data.table::fwrite(data.frame('common_genes' = common_genes),
                       file = paste0(out, '/model/', reference_name, '/common_genes_intersect.csv'))
    # filter each data set for common genes
    data.it = lapply(data, function(x){x[,common_genes]})
  } 
  else{
    common_genes = colnames(data[['ref']])
    # save common genes 
    data.table::fwrite(data.frame('common_genes' = common_genes),
                       file = paste0(out, '/model/', reference_name, '/common_genes_complete.csv'))
    # filter each data set for common genes
    data.it = lapply(data, function(x){
      ## Check which genes are missing
      miss.genes <- setdiff(common_genes,colnames(x))
      if(length(miss.genes) > 0){
        ## Create a matrix with zeros
        mtx.inp    <- matrix(0,
                             nrow = nrow(x),
                             ncol = length(miss.genes),
                             dimnames = list(rownames(x),
                                             miss.genes))
        x <- cbind(x,mtx.inp)
      }
      #Keep only the genes in the reference and reorder 
      x[,common_genes,drop=F]
    })
  }
  #----- SAVE DATA ----------------------------------------
  #Check if there is ANY cell with all zero values across features in reference
  if(any(rowSums(data.it[['ref']]) == 0)){
    stop("@ After processing the reference contains cells with zeros across all features. Remove them and re run the pipeline")
  }
  # save reference 
  
  tmp = data.it[['ref']] %>% rownames_to_column()
  colnames(tmp)[1] = " "
  

  data.table::fwrite(tmp,
                     file = paste0(out, '/model/', reference_name, '/expression_',ft_mth,'.csv'),
                     sep = ',')
  
  if(pipeline_mode == "annotation"){
  # save query 
  query_names = names(data.it)[!names(data.it) == 'ref']
  for(q in query_names){
    print(q)
    if(any(rowSums(data.it[[q]]) == 0)){
      stop(glue("@ After processing the query {q} contains cells with zeros across all features. Remove them and re run the pipeline"))
    } 
    tmp = data.it[[q]] %>% rownames_to_column()
    colnames(tmp)[1] = " "
    
    data.table::fwrite(tmp,
                       file = paste0(out, '/', q, '/', reference_name, '/expression_',ft_mth,'.csv'),
                       sep = ',')
  }
  }
  rm(data.it)
}
