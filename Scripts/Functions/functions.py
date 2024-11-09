
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
  import sys
    ## For each reference
  for ref in config['references'].keys():
      use_pretrain = False  # Default to False unless criteria are met
        # Check if pretrain models path is specified
      pretrain_models_path = config["references"][ref].get("pretrain_models")
      if pretrain_models_path:
          # Check if the path exists
          if not os.path.exists(pretrain_models_path):
              sys.exit("@ The pretrained models path doesn't exist.")
          # Verify tools with pretrained models
          tools_pretrain_models = [
              tool_pretrain for tool_pretrain in tools_run 
              if tool_pretrain in os.listdir(pretrain_models_path)
          ]
            # Set paths for expression matrix and labels
          config["references"][ref]["expression"] = os.path.join(pretrain_models_path, "expression_complete.csv")
          config["references"][ref]["labels"] = os.path.join(pretrain_models_path, "downsampled_labels.csv")
          config["references"][ref]['batch'] = None
          # Adjust reference parameters if necessary
          if config["references"][ref].get('convert_ref_mm_to_hg'):
              print("@ The original feature space of the models cannot be changed, setting convert_ref_mm_to_hg to False.")
              config["references"][ref]['convert_ref_mm_to_hg'] = False
          if int(config["references"][ref].get('min_cells_per_cluster')) > 0:
              print("@ The original feature space of the models cannot be changed, removing min_cells_per_cluster.")
              config["references"][ref]['min_cells_per_cluster'] = 0
          if config["references"][ref].get('downsample').get('value') > 0:
              print("@ The original feature space of the models cannot be changed, disabling downsampling.")
              config["references"][ref]['downsample']['value'] = 0
              config["references"][ref]['downsample']['stratified'] = False
            # Mark that pretrained models are being used
          use_pretrain = True
      # Save whether this reference uses pretrained models
      config["references"][ref]["use_pretrain"] = use_pretrain
      # Configure model directories and gene selection for each tool
      for tool in tools_run:
          config["references"][ref][tool] = {}
          if use_pretrain and (tool in tools_pretrain_models):
              # If using pretrained models, set the model directory and complete gene selection
              config["references"][ref][tool]["model_dir"] = pretrain_models_path
              config["references"][ref][tool]["gene_selection"] = "complete"
          else:
              # Otherwise, use the output directory and user-defined gene selection
              config["references"][ref][tool]["model_dir"] = os.path.join(config['output_dir'], "model", ref)
              if tool.startswith("scPred_"):
                  config["references"][ref][tool]["gene_selection"] = config.get("scPred", {}).get("gene_selection", "default")
              else:
                  config["references"][ref][tool]["gene_selection"] = config.get(tool, {}).get("gene_selection", "default")

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
