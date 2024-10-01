#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28

@author: tomas.vega
"""
## Python 3.11.2
#--------------- Libraries -------------------
import pandas as pd
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import anndata as ad
import sys
import os
import scanpy as sc
import numpy as np
# import pickle
import random

# Set seed
random.seed(123456) 

# Function to get the column name and maximum value for each row
def get_max_column_and_value(row):
    pred_label = row.idxmax()
    proba_label = row.max()
    return pred_label, proba_label
  
#--------------- Parameters -------------------
sample_path = str(sys.argv[1])
model_path = os.path.dirname(str(sys.argv[2]))
out_path = str(sys.argv[3])
out_other_path = os.path.dirname(str(sys.argv[3]))
threshold = float(sys.argv[4])

#--------------- Data -------------------

# read query matrix
print('@ READ QUERY')
query = pd.read_csv(sample_path,
                    index_col=0,
                    sep=',',
                    engine='c') ## could be pyarrow to make it faster but it needs to be install on the module and has many problems with dependencies

print('@ DONE')

# # load model 
# print('@ LOAD MODEL')
# scANVI_model = sca.models.SCANVI.load(dir_path= model_path)
# print('@ DONE')

# Query preprocessing
query = ad.AnnData(X = query,
                   obs = dict(obs_names=query.index.astype(str)),
                   var = dict(var_names=query.columns.astype(str))
                   )

## This step is necessary according to the authors
query = remove_sparsity(query)

query.obs['batch'] = 'query'
query.obs['labels'] = "Unknown"

#----------- Predict scANVI --------
## Need the intersection
model = sca.models.SCANVI.load_query_data(
    query,
    model_path,
    freeze_dropout = True,
)

model._unlabeled_indices = np.arange(query.n_obs)
model._labeled_indices = []

model.train(
    max_epochs=100,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)
# soft True to get the probability matrix (if I put false it only returns the label where it gets the highest prob)
pred_proba = model.predict(soft = True)                                   

# Create a new column 'max_column' with the column name containing the maximum value for each row
pred_proba['pred_label'], pred_proba['proba_label'] = zip(*pred_proba.apply(get_max_column_and_value, axis=1))

# Create a new column 'unknown_max_column' to store 'max_column' as 'unknown' if 'max_value' is lower than the threshold
pred_proba['pred_label_reject'] = pred_proba.apply(lambda row: 'Unknown' if row['proba_label'] < threshold else row['pred_label'], axis=1)


print('@ WRITTING PREDICTIONS')
pred_df = pd.DataFrame({'cell': pred_proba.index,
                       'scANVI': pred_proba.pred_label_reject})
pred_df.to_csv(out_path,
               index = False)
print('@ DONE')


#------------- Other outputs --------------

# Save the prob matrix
print('@ WRITTING PROB MATRIX ')
filename = out_other_path + '/scANVI_pred_score.csv'
pred_proba.to_csv(filename,
          index=True) #True because we want to conserve the rownames (cells)
print('@ DONE ')

# Save the prob matrix
print('@ SAVING THE QUERY MODEL')
filename = out_other_path + '/scANVI_model/'
model.save(filename, overwrite=True)
print('@ DONE ')
