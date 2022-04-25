#Parse Chad's reactivity data and output list of pairs 
import os
import sys
import pandas as pd
import numpy as np
from collections import defaultdict
import pickle
from itertools import product, chain


#Three sources of reactivity data:
#turnover for TropB, AzaH, SorbC
#Chirality data (stereo folder)
#Reactivity data for set of 8 ligands for first JGI set
#A model ligand pair is considered reactive if it shows reactivity in any of these sources

#Sources:
jgidata_file = "FMO_substrate_smiles_v2.xlsx"
turnover_file = 'turnover.py'

# pair dictionary: model -> set(reactive), set(all_ligands)
pair_dict = dict()

def add_pair(model, ligand, reactivity):
  if model not in pair_dict:
    pair_dict[model] = [set(), set()]
  
  #add ligand to set of all ligands
  pair_dict[model][1].add(ligand)
  
  #if reactive add to set of reactive ligands
  if reactivity:
    pair_dict[model][0].add(ligand)
  

#Turnover data:
from turnover import *

print('Turnover data:')
print('Ligand, TropB Turnover, AzaH Turnover, SorbC Turnover')
for lig, turnovers in turnover_data.items():
  print(lig, turnovers)

#add turnover data to pair_dict:
for ligand, turnovers in turnover_data.items():
  TropB_turnover, AzaH_turnover, SorbC_turnover = turnovers
  if TropB_turnover > 0:
    add_pair('tropb', ligand, True)
  else:
    add_pair('tropb', ligand, False)
    
  if AzaH_turnover > 0:
    add_pair('azah', ligand, True)
  else:
    add_pair('azah', ligand, False)
    
  if SorbC_turnover > 0:
    add_pair('sorbc', ligand, True)
  else:
    add_pair('sorbc', ligand, False)



#Chirality data:
#replace abbreivations with full name and remove Anc
chirality_ligands = ['18-2', '1-2', '17-2', 'ome', 'ketone']
for ligand in chirality_ligands:
  data_df = pd.read_excel(f'../../stereo/data/all_stereo_{ligand}.xlsx', header = 0, index_col = 0).fillna('')
  for model in data_df.index:
    stereo = data_df['R_frac'][model]
    if stereo:
      if stereo == '-':
        add_pair(model, ligand, False)
      else:
        add_pair(model, ligand, True)
  

#JGI data
jgi_df = pd.read_excel(jgidata_file, sheet_name = 1, index_col = 0, header=6).fillna('-')  #read actual reactivity
jgi_df.drop(jgi_df.filter(regex="Unname"),axis=1, inplace=True) #drop unnamed columns

print('JGI data:')
print(jgi_df)


for ligand in jgi_df.columns:
  for model in jgi_df.index:
    reactivity = jgi_df[ligand][model]
    model = str(model)
    ligand = str(ligand)
    
    if reactivity == '+':
      add_pair(model, ligand, True)
    else:
      add_pair(model, ligand, False)
      
      
#write model ligand pairs to txt files for train and test sets
pickle.dump(pair_dict, open('pair_dict.pkl', 'wb'))
pairs = list(chain.from_iterable( #flatten 
    [list(product([protein], ligands[1])) for protein, ligands in pair_dict.items()] #get protein ligand combinations
    )) 
pickle.dump(pairs, open('pairs.pkl', 'wb'))

print('Summary of reactivity data')

with open('pairs.txt', 'w') as pair_writer:
     
  models = list(pair_dict.keys())
  for model in models:
    reactive_ligands = pair_dict[model][0]
    unreactive_ligands = pair_dict[model][1].difference(pair_dict[model][0])
    if reactive_ligands:
      reactive_ligands = [x + '_+' for x in reactive_ligands]
      
    if unreactive_ligands:
      unreactive_ligands = [x + '_-' for x in unreactive_ligands]
      
    pair_writer.write(model + ' ' + ' '.join(reactive_ligands) + ' ' + ' '.join(unreactive_ligands) + '\n')
    print(model)
    print(reactive_ligands)
    print(unreactive_ligands)
      



