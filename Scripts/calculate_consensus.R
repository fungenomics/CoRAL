
library(tidyverse)
library(glue)


# get path to script to load custom functions  
initial.options = commandArgs(trailingOnly = FALSE)
file.arg.name = "--file="
script.name = sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)]) 

source(paste0(dirname(script.name), "/Functions/functions.R"))

set.seed(1234)

args = commandArgs(trailingOnly = T)
pred_path = args[1]
summary_path = args[2]
tools = strsplit(args[3], split = ' ')[[1]]
consensus_tools = strsplit(args[4], split = ' ')[[1]]
consensus_type = args[5]
min_agree = as.numeric(strsplit(args[6], split = ' ')[[1]])
ontology_column = args[7]
ontology_path = args[8]
metrics_path = args[9]
CAWPE_type = strsplit(args[10], split = ' ')[[1]]
alpha = as.numeric(strsplit(args[11], split = ' ')[[1]])
metric = args[12]

print(CAWPE_type)


if(consensus_tools[1] == 'all'){
  consensus_tools = tools
}

print('CONSENSUS TOOLS: ')
print(consensus_tools)

print('CONSENSUS TYPE: ')
print(consensus_type)

# read predictions 
consensus = read_prediction_files(pred_path, tools = tools)

print(head(consensus))

# read ontology
ontology = data.table::fread(ontology_path,
                             sep = ",",
                             header = T) %>% as.data.frame

check_ontology_hierarchy(ontology)

print(head(ontology))

# Majority Vote 
if(consensus_type == 'majority' & all(min_agree != 0)){
  consensus[,tools] <- lapply(consensus[, tools, drop=F],
                          FUN = function(x){
                            apply_ontology(df_ontology = ontology,
                                           pred = x,
                                           from = 'label',
                                           to = ontology_column)
                          })

  # harmonize unresolved values from the different tools so that they are the same 
  tmp = harmonize_unresolved(consensus, unique(ontology[[ontology_column]])) %>% 
        select(all_of(consensus_tools))

  for(mn_ag in min_agree){
    consensus[,paste0("Consensus_", as.character(mn_ag))] <- apply(tmp, 1, get_consensus, min = mn_ag)
  }

  # Add Max Vote and Entropy to table (QC)
  consensus[,"majority_MaxVote"] = apply(tmp, 1, get_max)
  consensus[,"majority_Entropy"] = apply(tmp, 1, get_entropy)

  rm(tmp)

}else if(consensus_type == "CAWPE" & all(alpha != 0)){

   # get the paths to the probability files 
   prob_files = list.files(pred_path, 
                          pattern = paste0(paste0(consensus_tools, '_pred_score.csv'), collapse = '|'), 
                          recursive = T, 
                          full.names = T)
    
   for(CW_tp in CAWPE_type){
    # get the column names to group by (just tool or tool + label)
    cols = get_cawpe_columns(CW_tp)
    
    # read metrics table
    print(metric)
    metrics_file = data.table::fread(metrics_path,
                                     check.names = F,
                                     sep = ",")

    print(head(metrics_file))
    
    # check if the metric sepcified is available 
    if(!metric %in% colnames(metrics_file)){
       avail_met = paste0(colnames(metrics_file)[!colnames(metrics_file) %in% c('tool', 'fold')], collapse=', ')
       stop(paste0("@ metric ", metric, " not in file. Choose one of the following: ", avail_met))
    } 

    print(cols)
    print(consensus_tools)

    # check if the tools sepcified have the metrics calculated
    if(!all(consensus_tools %in% unique(metrics_file$tool))){
      ms_tool = paste0(consensus_tools[!(consensus_tools %in% unique(metrics_file$tool))], collapse=', ')
      stop(paste0("@",metric, " Metric from Tools ", ms_tool, " is missing"))
    }
    # read metric file and caluclate mean metric per tool or tool + label 
    metrics_file = metrics_file %>% 
      filter(tool %in% consensus_tools) %>%
      mutate(ontology = class) %>% 
      mutate(metric = .data[[metric]]) %>% 
      group_by(across(all_of(cols))) %>%
      summarize(mean_metric = mean(metric)) %>%
      ungroup()
     
    print(head(metrics_file))
    
    for(aph in alpha){
      # read probability matrixes
      data = lapply(prob_files, function(f){data.table::fread(f, check.names = F) %>% 
            rename(cellname = V1) %>% 
            mutate(tool = basename(dirname(f)))})
      
      # merge prob matrices and join with metric scores
      data = bind_rows(data) %>% 
            select(cellname, all_of(unique(ontology[['label']])), tool) %>%
            pivot_longer(!c('cellname', 'tool'), 
                        names_to = 'class', 
                        values_to = 'prob') %>% 
        as.data.frame()
      
      data$ontology = apply_ontology(df_ontology = ontology,
                                     pred = data$class,
                                     from = 'label',
                                     to = ontology_column)  
      # I treat everything like an ontology, for the prob I sum the prob,
      # for the metric since it is already calculated at the ontology level, 
      # just do unique
      data <- data  %>%
            left_join(metrics_file, by = cols) %>%
            na.omit() %>%
            group_by(cellname,tool, ontology) %>% 
            summarise(prob_ont = sum(prob),
                      mean_metric = unique(mean_metric))
      # calculate cawpe       
      data$CAWPE = (as.numeric(data$mean_metric)^aph) * as.numeric(data$prob_ont)
         
      # add up the CAWPE scores for each class for each cell 
      data = data %>%
           group_by(cellname, ontology) %>%
        #Doing the mean it norm by the number of tools and is comparable between ref (value between 0 and 1)
           summarise(CAWPE = mean(CAWPE)) %>% 
           ungroup() 
      
      data = data %>% 
             group_by(cellname) %>% 
             mutate(prop = CAWPE / sum(CAWPE)) %>% 
             mutate(entropy = ifelse(test = length(ontology) > 1,
                                     yes = (-sum(prop * log2(prop),na.rm = T))/log2(length(ontology)),
                                     no = 0)) %>% 
            ungroup() %>% select(-prop) %>% 
            as.data.frame()
        
      data.wd <- data %>% 
                 pivot_wider(names_from = ontology,
                             values_from=CAWPE,
                             values_fill = 0)
        
      # save CAWPE scores 
      data.table::fwrite(data.wd, 
                         file = paste0(dirname(summary_path),'/',CW_tp,'_', aph, '_', ontology_column, '_scores.csv' ))
      
      rm(data.wd)
      
      data = data %>% 
             group_by(cellname) %>%
             dplyr::slice(which.max(CAWPE)) %>%
             rename(Consensus = ontology)  %>% 
             as.data.frame()
  
      rownames(data) = data$cellname
  
      print(head(data))
        
      consensus[, paste0("CAWPE_",CW_tp,"_",aph)]         = data[consensus$cellname, "CAWPE", drop=T]
      consensus[, paste0("CAWPE_entropy_",CW_tp,"_",aph)] = data[consensus$cellname, "entropy", drop=T]
      consensus[, paste0("Consensus_",CW_tp,"_",aph)]     = data[consensus$cellname, "Consensus", drop=T]

     
     }
   }
}

print("@Modifying the original tool labels to ontology")
consensus[,tools] = lapply(consensus[,tools,drop=F],
                           FUN = function(x){
                             apply_ontology(df_ontology = ontology,
                                            pred = x,
                                            from = 'label',
                                            to = ontology_column)
                           }
)

data.table::fwrite(consensus, summary_path, sep = '\t')

#--------------------------------------------------
