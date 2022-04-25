#For a list of proteins and ligands get a list of structures and docked poses
import os
import sys
import shutil
import pandas as pd
import pickle

proteins = [
    'tropb',
    'afod',
    'azah',
    '312',  #TROPB ANCESTOR
    '321',  #TROPB ANCESTOR
    '331',  #TROPB ANCESTOR
    '332',  #TROPB ANCESTOR
    '369',  #AFOD ANCESTOR
    '369a', #AFOD ANCESTOR
    '370',  #AFOD ANCESTOR
    '371',  #AFOD ANCESTOR
    '373',  #AFOD ANCESTOR
    '374',  #AFOD ANCESTOR
]

ligands = ['1-2',  #3MO 
           '18-2', #MODEL SUBSTRATE
           ]

#reactivity pairs:
pairs = pickle.load(open('../../anc_model/reactivity/data/pairs.pkl', 'rb'))
proteins = set()
ligands = set()
for protein, ligand in pairs:
    proteins.add(protein)
    ligands.add(ligand)

study = 'prot_1_fftdock_3' #fftdock poses generated in grid with eps 3 and minimized in explicit protein in eps 1


#Get all alphafold structures of enzymes:
structure_path = '/home/azamh/anc_model_old/add_cofactor/pdb_with_fad'
print('GETTING ALPHAFOLD STRUCTURES')

for protein in proteins:
    protein_with_fad = os.path.join(structure_path, f'{protein}_fad.pdb')
    if not os.path.isfile(protein_with_fad):
        print(f'\tNO STRUCTURE FOR {protein}')
        continue
    else:
        shutil.copy(protein_with_fad, './pdb_with_fad')
        print(f'\tSTRUCTURE FOR {protein} COPIED')
        
#Pull all docked poses
pose_path = '/home/azamh/anc_model/docking/pdb_top_15/{protein}_{ligand}/{protein}_{ligand}_{study}_0.pdb'
print('\nGETTING LIGAND DOCKED POSES')

for protein in proteins:
    for ligand in ligands:
        pose = pose_path.format(protein = protein, ligand = ligand, study = study)
        
        if not os.path.isfile(pose):
            print(f'\tNO DOCKED POSE FOR {protein} and {ligand}')
        else:
            shutil.copy(pose, './top_poses')
            print(f'\tDOCKED POSE FOR {protein} and {ligand} COPIED')
        
#Print true and predicted stereochemistries from docking:
pred_df = pd.read_excel('/home/azamh/anc_model/stereo/pred_stereo/pred_stereo.xlsx', header = 0, index_col = 0)
pred_df = pred_df.fillna('')
cluster_path = '/home/azamh/anc_model/docking/cluster/{protein}_{ligand}_{study}.xlsx'

copy_df = pd.DataFrame(columns = ['PROTEIN', 'LIGAND', 'TRUE', 'PRED', "ENERGY"])
print('\nPRINTING PREDICTED AND TRUE STEREOCHEMISTRIES (R FRACTION)')

for protein in proteins:
    for ligand in ligands:
        key = f"('{protein}', '{ligand}')"
        if protein == 'afod_crystal':
            key = f"('afod', '{ligand}')"
        if key in pred_df.index:
            true_stereo = pred_df['True'].loc[key]
            pred_stereo = pred_df[study].loc[key]
            
            if true_stereo or pred_stereo:
                
                #add in docked pose energy from clusters
                dock_ener = 'NO POSE'
                cluster_file = cluster_path.format(protein = protein, ligand = ligand, study = study)
                if os.path.isfile(cluster_file):
                    cluster_df = pd.read_excel(cluster_file, header=0, index_col=0)
                    dock_ener = cluster_df['min_ener'][0]
                
                
                data = {"PROTEIN":protein,
                        "LIGAND":ligand,
                        "TRUE":true_stereo,
                        "PRED":pred_stereo,
                        "ENERGY":dock_ener,
                }
                copy_df = copy_df.append(data, ignore_index=True)


                
print(copy_df)
copy_df.to_excel('stereo_energy.xlsx')


sys.exit()
        
