library(tidyverse)

set.seed(1234)

args = commandArgs(trailingOnly = TRUE)

out = args[1]
lab_path = args[2]
ontology_path = args[3]                                                                                                                            
ontology_columns = strsplit(args[4], split = ' ')[[1]]

print(out)
print(lab_path)
print(ontology_path)
print(ontology_columns)

#----- SAVE ONTOLOGY -----------------------------------

dir.create(dirname(out),
           recursive = T)

lab = data.table::fread(lab_path,
                        sep = ",",
                        header = T) %>% as.data.frame()

print(lab)

if(length(ontology_columns) == 1 & ontology_columns[1] == 'label'){
  
  ont = data.frame(label = unique(lab$label))
  print(ont)

  data.table::fwrite(ont,
                     file = out,
                     sep = ',')
}else{
  ont = data.table::fread(ontology_path,
                          sep = ",",
                          header = T) %>% as.data.frame() 
  
  #Filter from the ontology the labels that were removed from the reference in the preprocess 
  ont <- ont[ont$label %in% lab$label,,drop=F]
  print(ont)
  data.table::fwrite(ont,
                     file = out,
                     sep = ',')
}

#--------------------------------------------------------
