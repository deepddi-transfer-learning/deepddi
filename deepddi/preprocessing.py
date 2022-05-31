import os
import glob
import pickle
import copy
import argparse
import itertools
import pandas as pd
import time
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit import Chem
import tqdm

def parse_food_input(input_file):
    parsed=open('data/parsed_input.txt','w+')
    parsed.write('Prescription	Drug name	Smiles\n')
    food_compound = pd.read_csv('database/food_compounds_lookup.csv')
    merged=pd.read_csv('database/drug_info_combined.csv')
    all_d=[]
    all_f=[]

    first_line=True
    with open(input_file, 'r') as fp:
        for line in fp:
            each_input=[i.lower() for i in line.strip().split('\t')]
            if first_line:
                for i in each_input:
                    all_d.append(i)
                first_line=False                      
            else:
                for i in each_input:
                    all_f.append(i)
    
    assert len(all_d)>=1 and len(all_f)>=1, 'No valid pairs entered'
    count=0
    for i in all_d:
        drug_name,drug_com=i.strip().lower().split('|')
        name_i = drug_name+'('+drug_com+')'+'\t'
        find_drug_1=merged.loc[merged['Name'].str.lower()==drug_com.lower()]
        smile_i=find_drug_1['Smiles'].values[0]+'\n'

        for j in all_f:
            food_name,food_comp=j.strip().lower().split('|')
            
            each_i=str(count)+'\t'
            each_j=str(count)+'\t'
        #     find the food in the approved food list
            
            each_i += name_i
            each_j += food_name+'('+food_comp+')'+'\t'
            
            find_food_2 = food_compound.loc[food_compound['name'].str.lower()==food_comp.lower()]
            smile_j = find_food_2['cas_number'].values[0]+'\n'
            
            each_i+=smile_i
            each_j+=smile_j
            parsed.write(each_i)
            parsed.write(each_j)
            count+=1
    parsed.close()
    return 


    

def parse_drug_input(input_file):
    parsed=open('data/parsed_input.txt','w+')
    parsed.write('Prescription	Drug name	Smiles\n')
    merged=pd.read_csv('data/Drug_info_combined.csv')
    approved_drugs=set(merged['Name'].str.lower())
    current_drug=[]
    other_drugs=[]

    first_line=True
    with open(input_file, 'r') as fp:
        for line in fp:
            if first_line:
                current_drug=line.strip().split('\t')
                first_line=False
            else:
                other_drugs=line.strip().split('\t')

    count=0
    for i in current_drug:
        drug_name_i,drug_com_i=i.strip().lower().split('|')
        name_i = drug_name_i+'('+drug_com_i+')'+'\t'
        find_drug_1=merged.loc[merged['Name'].str.lower()==drug_com_i]

        for j in other_drugs:
            drug_name_j,drug_com_j=j.strip().lower().split('|')

            each_i=str(count)+'\t'
            each_j=str(count)+'\t'
        #     find the drug in the approved drug list
            find_drug_2=merged.loc[merged['Name'].str.lower()==drug_com_j]

            name_j = drug_name_j+'('+drug_com_j+')'+'\t'

            each_i += name_i
            each_j += name_j
            
            smile_i=find_drug_1['Smiles'].values[0]+'\n'
            smile_j=find_drug_2['Smiles'].values[0]+'\n'
            
            each_i+=smile_i
            each_j+=smile_j
            parsed.write(each_i)
            parsed.write(each_j)
            count+=1
    parsed.close()
    return 
            

def parse_DDI_input_file(input_file, output_file):
    drug_pair_info = {}
    drug_smiles_info = {}
    with open(input_file, 'r') as fp:
        for line in fp:
            sptlist = line.strip().split('\t')
            prescription = sptlist[0].strip()
            drug_name = sptlist[1].strip()
            smiles = sptlist[2].strip()
            
            drug_smiles_info[(prescription, drug_name)] = smiles
            if prescription not in drug_pair_info:
                drug_pair_info[prescription] = [drug_name]
            else:
                drug_pair_info[prescription].append(drug_name)
            
        
    out_fp = open(output_file, 'w')
    for each_prescription in drug_pair_info:
        drug_names = drug_pair_info[each_prescription]
        for each_set in itertools.combinations(drug_names, 2):
            drug1 = each_set[0].strip()
            drug1_smiles = drug_smiles_info[(each_prescription, drug1)] 
            
            drug2 = each_set[1].strip()
            drug2_smiles = drug_smiles_info[(each_prescription, drug2)] 
            out_fp.write('%s\t%s\t%s\t%s\t%s\n'%(each_prescription, drug1, drug1_smiles, drug2, drug2_smiles))
    out_fp.close()
    return

def calculate_drug_similarity(drug_dir, input_dir, output_file):
    drugbank_drugs = glob.glob(drug_dir + '*')
    input_drugs = glob.glob(input_dir + '*')
    drug_similarity_info = {}
    for each_drug_id1 in drugbank_drugs:
        drugbank_id = os.path.basename(each_drug_id1).split('.')[0]
        drug_similarity_info[drugbank_id] = {}
        drug1_mol = Chem.MolFromMolFile(each_drug_id1)
        drug1_mol = AllChem.AddHs(drug1_mol)
        for each_drug_id2 in input_drugs:            
            input_drug_id = os.path.basename(each_drug_id2).split('.')[0]
            drug2_mol = Chem.MolFromMolFile(each_drug_id2)
            drug2_mol = AllChem.AddHs(drug2_mol)
            fps = AllChem.GetMorganFingerprint(drug1_mol, 2)
            fps2 = AllChem.GetMorganFingerprint(drug2_mol, 2)
            score = DataStructs.TanimotoSimilarity(fps, fps2)
            drug_similarity_info[drugbank_id][input_drug_id] = score

    df = pd.DataFrame.from_dict(drug_similarity_info)
    df.to_csv(output_file)

    
def calculate_structure_similarity(drug_dir, input_file, output_file, drug_list):
    drugbank_drugs = glob.glob(drug_dir + '*')
    all_input_drug_info = {}
    with open(input_file, 'r')as fp:
        for line in fp:
            sptlist = line.strip().split('\t')
            prescription = sptlist[0].strip()
            drug1 = sptlist[1].strip()
            smiles1 = sptlist[2].strip()
            drug2 = sptlist[3].strip()
            smiles2 = sptlist[4].strip()
            if drug1 not in all_input_drug_info:
                all_input_drug_info[drug1] = smiles1
            if drug2 not in all_input_drug_info:
                all_input_drug_info[drug2] = smiles2            
        
    drug_similarity_info = {}
    for input_drug_id in all_input_drug_info:   
        try:
            each_smiles = all_input_drug_info[input_drug_id]
            drug2_mol = Chem.MolFromSmiles(each_smiles)
            drug2_mol = AllChem.AddHs(drug2_mol)            
        except:
            continue
        drug_similarity_info[input_drug_id] = {}
        for each_drug_id1 in drugbank_drugs:            
            drugbank_id = os.path.basename(each_drug_id1).split('.')[0]
            
            drug1_mol = Chem.MolFromMolFile(each_drug_id1)        
            drug1_mol = AllChem.AddHs(drug1_mol)    
            
            fps = AllChem.GetMorganFingerprint(drug1_mol, 2)
            fps2 = AllChem.GetMorganFingerprint(drug2_mol, 2)
            score = DataStructs.TanimotoSimilarity(fps, fps2)
            drug_similarity_info[input_drug_id][drugbank_id] = score
            
    df = pd.DataFrame.from_dict(drug_similarity_info)
    df = df.T
    df = df[drug_list]
    df.to_csv(output_file)

def calculate_pca(similarity_profile_file, output_file, pca_model):
    with open(pca_model, 'rb') as fid:
        pca = pickle.load(fid)
        df = pd.read_csv(similarity_profile_file, index_col=0)

        X = df.values
        X = pca.transform(X)

        new_df = pd.DataFrame(X, columns=['PC_%d' % (i + 1) for i in range(50)], index=df.index)
        new_df.to_csv(output_file)

def generate_input_profile(input_file, pca_profile_file):    
    df = pd.read_csv(pca_profile_file, index_col=0)
    # df.index = df.index.map(str)
    
    all_drugs = []
    interaction_list = []
    with open(input_file, 'r') as fp:
        for line in fp:
            sptlist = line.strip().split('\t')
            prescription = sptlist[0].strip()
            drug1 = sptlist[1].strip()
            drug2 = sptlist[3].strip()
            all_drugs.append(drug1)
            all_drugs.append(drug2)
            if drug1 in df.index and drug2 in df.index:
                interaction_list.append([prescription, drug1, drug2])
                interaction_list.append([prescription, drug2, drug1])
    
    drug_feature_info = {}
    columns = ['PC_%d' % (i + 1) for i in range(50)]
    for row in df.itertuples():
        drug = row.Index
        feature = []
        drug_feature_info[drug] = {}
        for col in columns:
            val = getattr(row, col)
            feature.append(val)
            drug_feature_info[drug][col] = val

    new_col1 = ['1_%s'%(i) for i in columns]
    new_col2 = ['2_%s'%(i) for i in columns]
    
    DDI_input = {}
    for each_drug_pair in tqdm.tqdm(interaction_list):
        prescription = each_drug_pair[0]
        drug1 = each_drug_pair[1]
        drug2 = each_drug_pair[2]
        key = '%s_%s_%s' % (prescription, drug1, drug2)
        
        DDI_input[key] = {}
        for col in columns:
            new_col = '1_%s'%(col)
            DDI_input[key][new_col] = drug_feature_info[drug1][col]
            
        for col in columns:
            new_col = '2_%s'%(col)
            DDI_input[key][new_col] = drug_feature_info[drug2][col]

    new_columns = []
    for i in [1,2]:
        for j in range(1, 51):
            new_key = '%s_PC_%s'%(i, j)
            new_columns.append(new_key)
            
    df = pd.DataFrame.from_dict(DDI_input)
    df = df.T
    df = df[new_columns]
    # df.to_csv(output_file)
    return df