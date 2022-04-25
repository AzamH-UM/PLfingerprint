#Look at distribution energy vs reactivity
import os
import pandas as pd
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
import math
import numpy as np

pairs = pickle.load(open('../pairs.pkl', 'rb'))
pair_dict = pickle.load(open('../pair_dict.pkl', 'rb'))  #protein -> reactive ligands, all ligands
ligands = list({pair[1] for pair in pairs})
train_ligs = ['1-2', '2-2', '17-2', '18-2', '32-2', 'a1', 'a2', 'ome', 'ketone']
study = 'prot_1_fftdock_3'

data_df = pd.DataFrame(columns = ['Protein', 'Ligand', 'Reactivity', 'Energy'])

for protein, ligand in pairs:
    if ligand in train_ligs:
        clusterfile = f'../../../docking/cluster/{protein}_{ligand}_{study}.xlsx'
        if os.path.isfile(clusterfile):
            print(protein, ligand)
            energy = np.nan
            reactivity = np.nan
            cluster_df = pd.read_excel(clusterfile, header = 0, index_col = 0)
            energy = cluster_df['min_ener'][0]
            
            if energy > 100000 or energy < -100000:
                continue
            
            reactive_ligands, all_ligands = pair_dict[protein]
            if ligand in reactive_ligands:
                reactivity = 'reactive'
            elif ligand in all_ligands:
                reactivity = 'unreactive'
            else:
                raise ValueError()
            
            
            data = {'Protein':protein,
                    'Ligand':ligand,
                    'Reactivity':reactivity,
                    'Energy':energy}
            
            data_df = data_df.append(data, ignore_index=True)
    

data_df.to_excel('data.xlsx')
plt.figure(figsize= (10,5))
plt.xlabel('Energy')
plt.ylabel('Ligand')
plt.ylim(-100, 50)
sns.stripplot(data= data_df, x = 'Ligand', y = 'Energy', hue = 'Reactivity')
plt.legend()
plt.savefig('stripplot.png')