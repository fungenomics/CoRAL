
# set parameters related to converting gene names 
def set_gene_conversion_parameters(config):
  for ref in config['references'].keys():
    try:
      config["references"][ref]['convert_ref_mm_to_hg']
    except:
      config["references"][ref]['convert_ref_mm_to_hg'] = False

# set tools-ref parameters
# the parameters of the output_dir for atool and the gene_filtering 
# is specific for each reference, this allows to make more flexible
# and be able to use pretrain models.
# this could be extended to the other parameters if needed
def set_reference_pretrained_parameter(config, tools_run):
  import os
  ## For each reference
  for ref in config['references'].keys():
    # Check if the pretrain models were required
    try: 
      config["references"][ref]["pretrain_models"]
      # If so
      #Check if the path exist
      if not os.path.exists(config["references"][ref]["pretrain_models"]):
        sys.exit("@ The pretrained models path doesn't exist.")
      ## Then for the tools required to run, check if they have pretrain models
      tools_pretrain_models = [tool_pretrain for tool_pretrain in tools_run if tool_pretrain in os.listdir(config["references"][ref]["pretrain_models"])]
      # and flag this reference as a reference that is going to use pretrained models
      use_pretrain = True
    except:
      use_pretrain = False  
    finally:
      config["references"][ref]["use_pretrain"] = use_pretrain
      #For all the tools to run
      for tool in tools_run:
        # Create the tool-reference 
        config["references"][ref][tool] = {}
        # If the tool has a pretrain model the input use for the model is the pretrain path
        # and should use the complete method for the queries
        if use_pretrain and (tool in tools_pretrain_models):
        # if True:
          config["references"][ref][tool]["model_dir"] = config["references"][ref]["pretrain_models"]
          # if tool.startswith("scPred_"):
          #   config["references"][ref]["scPred"]["gene_selection"] = "complete"
          # else:
          config["references"][ref][tool]["gene_selection"] = "complete"
        #if don't the input is the one set by the pipeline and the gene_selection is the one specified by the user
        else:
          config["references"][ref][tool]["model_dir"] = config['output_dir'] + "/model/" + ref
          if tool.startswith("scPred_"):
            config["references"][ref][tool]["gene_selection"] = config["scPred"]["gene_selection"]
          else:
            config["references"][ref][tool]["gene_selection"] = config[tool]["gene_selection"]
      
# set parameters related to batch 
def set_reference_batch_parameters(config):
  for ref in config['references'].keys():
    try:
      config["references"][ref]['batch']
    except:
      config["references"][ref]['batch'] = None

# set parameters related to benchmarking directory 
def set_benchmark_directory(config,mode):
  import sys 
  for ref in config['references'].keys():
    try:
      config["references"][ref]['output_dir_benchmark']
    except:
      if mode == 'annotation':
        if (config["consensus"]["type"]["CAWPE"]["alpha"][0] != 0) & (config["consensus"]["type"]["CAWPE"]["mode"] != ""):
          sys.exit("@ For running CAWPE the output directory of the benchmarking pipeline should be specified")
        else:
          config["references"][ref]['output_dir_benchmark'] = ""
      else:
        sys.exit("@ In the benchmarking pipeline, all the directories should be specified")
# set parameters related to downsampling 
def set_downsampling_parameters(config):
  for ref in config['references'].keys():
    try: 
      config["references"][ref]['min_cells_per_cluster']
    except: 
      # 50 is the min cell per cluster by default
      config["references"][ref]['min_cells_per_cluster'] = 50

    try: 
      config["references"][ref]['downsample']['value']
    except: 
      config["references"][ref]['downsample'] = {}
      config["references"][ref]['downsample']['value'] = 0

    try:
      config["references"][ref]['downsample']['stratified']
    except:	
      config["references"][ref]['downsample']['stratified'] = False


# set parameters related to ontology 
def set_ontology_parameters(config,mode):
  import pandas as pd
  
  for ref in config['references'].keys():
    
    try:
      config["references"][ref]["ontology"]["ontology_path"]
      
      try:
        config["references"][ref]["ontology"]["ontology_column"]

      except:
        columns = pd.read_csv(config["references"][ref]["ontology"]["ontology_path"], nrows=0).columns.tolist()
        config["references"][ref]["ontology"]["ontology_column"] = columns[1:]

    except: 
      config["references"][ref]["ontology"] = {}
      if mode == 'annotation':
        config["references"][ref]["ontology"]["ontology_path"] = config["output_dir"] + "/model/" + ref + "/ontology/ontology.csv"
      else:
        config["references"][ref]["ontology"]["ontology_path"] = config["references"][ref]['output_dir_benchmark'] + "/" + ref + '/ontology/ontology.csv'
      config["references"][ref]["ontology"]["ontology_column"] = ['label']
      continue 
    
    if not isinstance(config["references"][ref]["ontology"]["ontology_column"], list): 
      config["references"][ref]["ontology"]["ontology_column"] = [config["references"][ref]["ontology"]["ontology_column"]]
    
    config["references"][ref]["ontology"]["ontology_column"] = ['label'] + config["references"][ref]["ontology"]["ontology_column"]

# return consensus methods 
def get_consensus_methods(config,mode= "None"): 
  import sys 
  if mode == 'benchmarking':
    if config["consensus"]["type"]["majority"]["min_agree"][0] == 0:
      sys.exit("@ The minimum number of tools required to agree in the majority consensus should be specified")
  else:
    consensus_run = []
    if config["consensus"]["type"]["majority"]["min_agree"][0] != 0:
      consensus_run.append("majority")
    
    if (config["consensus"]["type"]["CAWPE"]["alpha"][0] != 0) & (config["consensus"]["type"]["CAWPE"]["mode"] != ""):
      consensus_run.append("CAWPE")
    
    if len(consensus_run) == 0:
      sys.exit("@ At least one consensus type should be specified")
  
    return consensus_run

# return the tools to run adding the models to scPred
def get_tools_to_run(config, mode = "None"):
  tools_run = config['tools_to_run']
  if "scPred" in tools_run:
    method = config['scPred']['classifier']
    if not isinstance(method, list):
        #Convert into a list
        method = [method]
    tools_to_run = [tool + "_" + m if tool == "scPred" else tool for tool in tools_run for m in (method if tool == "scPred" else [""])]
  else:
    tools_to_run = tools_run
  if mode == "pretrain" and any(tool in tools_to_run for tool in ["scNym","scAnnotate","scID"]):
    import sys
    sys.exit("scNym/scAnnotate/scID cannot be run in the pretrain pipeline (it cannot be split into train and test). Please remove it from your config file and rerun the pipeline.")
  return(tools_to_run)

# return the tools to run adding the models to scPred
def get_consensus_tools(config):
  consensus_to_run = config['consensus']['tools']
  if 'all' != consensus_to_run:
    if "scPred" in consensus_to_run:
      method = config['scPred']['classifier']
      if not isinstance(method, list):
          #Convert into a list
          method = [method]
      consensus_to_run = [tool + "_" + m if tool == "scPred" else tool for tool in consensus_to_run for m in (method if tool == "scPred" else [""])]
  return(consensus_to_run)

# return the feature selection methods specified
def get_feature_selection_method(config,tools_run ):
  feature_selection_methods = []
  for ref in config['references'].keys():
    for tool in tools_run:
      if config["references"][ref][tool]['gene_selection'] not in feature_selection_methods:
        feature_selection_methods.append(config["references"][ref][tool]['gene_selection'])
  return(feature_selection_methods)
