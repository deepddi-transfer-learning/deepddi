import os
import glob
import copy
import argparse
import time
import pandas as pd

from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
import tqdm
pd.set_option('mode.chained_assignment', None)


def read_information_file(information_file):    
    interaction_info = {}
    fp = open(information_file, 'r')
    fp.readline()
    for line in fp:
        sptlist = line.strip().split('\t')
        interaction_type = sptlist[0].strip()
        sentence = sptlist[1].strip()
        interaction_info[interaction_type] =  sentence
    fp.close()
    return interaction_info


def read_drug_information(drug_information_file):
    drug_information = {}
    with open(drug_information_file, 'r') as fp:
        for line in fp:
            sptlist = line.strip().split('\t')
            drugbank_id = sptlist[0].strip()
            drugbank_name = sptlist[1].strip()
            action = sptlist[7].strip()
            pharmacological_action = sptlist[8].strip()
            target = sptlist[5].strip()
            
            if action != 'None' and target != 'None':                
                if drugbank_id not in drug_information:
                    drug_information[drugbank_id] = [target]
                else:
                    drug_information[drugbank_id].append(target)
    return drug_information

def read_drug_enzyme_information(drug_enzyme_information_file):
    drug_information = {}
    with open(drug_enzyme_information_file, 'r') as fp:
        for line in fp:
            sptlist = line.strip().split('\t')
            drugbank_id = sptlist[0].strip()
            uniprot_id = sptlist[4].strip()
            action = sptlist[5].strip()
            
            if uniprot_id != 'None' and action != 'None':
                if drugbank_id not in drug_information:
                    drug_information[drugbank_id] = [uniprot_id]
                else:
                    drug_information[drugbank_id].append(uniprot_id)
                    
    return drug_information

def read_known_DDI_information(known_DDI_file):
    left_ddi_info = {}
    right_ddi_info = {}
    with open(known_DDI_file, 'r') as fp:
        fp.readline()
        for line in fp:
            sptlist = line.strip().split('\t')

            left_drug = sptlist[0].strip()
            right_drug = sptlist[1].strip()
            interaction_type = sptlist[2].strip()
            
            if interaction_type not in left_ddi_info:
                left_ddi_info[interaction_type] = [left_drug]
            else:
                left_ddi_info[interaction_type].append(left_drug)
                
            if interaction_type not in right_ddi_info:
                right_ddi_info[interaction_type] = [right_drug]
            else:
                right_ddi_info[interaction_type].append(right_drug)            
    
    for each_interaction_type in left_ddi_info:
        left_ddi_info[each_interaction_type] = list(set(left_ddi_info[each_interaction_type]))
    
    for each_interaction_type in right_ddi_info:
        right_ddi_info[each_interaction_type] = list(set(right_ddi_info[each_interaction_type]))
        
    return left_ddi_info, right_ddi_info

def read_similarity_file(similarity_file):
    similarity_df = pd.read_csv(similarity_file, index_col=0)

    return similarity_df

def get_side_effects(df, target_drug, frequency=10):
    new_df = df[df['Drug name']==target_drug]
    new_df = new_df[new_df['MEAN'] >= frequency]
    drug_side_effect_info = {}
    
    for each_drug, each_df in new_df.groupby('Drug name'):
        side_effects = each_df['SIDE EFFECT']
        string_list = []
        for each_index, each_df in each_df.iterrows():
            side_effect = each_df['SIDE EFFECT']
            mean_frequency = each_df['MEAN']
            string_list.append('%s(%.1f%%)'%(side_effect, mean_frequency))
            
        drug_side_effect_info[each_drug] = ';'.join(string_list)
        
    string_list = []
    for each_drug in drug_side_effect_info:
        string_list.append('%s'%(drug_side_effect_info[each_drug]))
        
    side_effect_string = ';'.join(string_list)
    return side_effect_string

def read_side_effect_info(df, frequency=10):
    new_df = df[df['MEAN'] >= frequency]
    drug_side_effect_info = {}
    for each_drug, each_df in new_df.groupby('Drug name'):
        drug_side_effect_info[each_drug] = None

        string_list = []
        for each_index, each_df in each_df.iterrows():
            side_effect = each_df['SIDE EFFECT']
            mean_frequency = each_df['MEAN']
            string_list.append('%s(%.1f%%)'%(side_effect, mean_frequency))
            
        drug_side_effect_info[each_drug] = ';'.join(string_list)
    
    return drug_side_effect_info

def annotate_DDI_results(DDI_output_file, drug_information_file, drug_enzyme_information_file, similarity_file, 
known_DDI_file, output_file, side_effect_information_file, model_threshold, structure_threshold):
    drug_information = read_drug_information(drug_information_file)    
    drug_enzyme_information = read_drug_enzyme_information(drug_enzyme_information_file)    
    
    left_ddi_info, right_ddi_info = read_known_DDI_information(known_DDI_file)    
    similarity_df = read_similarity_file(similarity_file)
    DDI_prediction_df = pd.read_csv(DDI_output_file, sep='\t')
    side_effect_df = pd.read_csv(side_effect_information_file, sep='\t')
    drug_side_effect_info = read_side_effect_info(side_effect_df, frequency=10)
    fp = open(output_file, 'w')
    fp.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%('Prescription','Drug_pair', 'Interaction_type',
     'Sentence', 'DDI_prob', 'DDI_prob_std', 'Confidence_DDI',   
     'Side effects (left)', 'Side effects (right)', 'Similar approved drugs (left)', 'Similar approved drugs (right)',
     'drug1', 'drug2'))

    s = time.time()
    for row in tqdm.tqdm(DDI_prediction_df.itertuples(), total=len(DDI_prediction_df)):

        prescription = str(getattr(row, 'Prescription'))
        drug_pair = getattr(row, 'Drug_pair')
        
        left_drug, right_drug = drug_pair.split('_')
        DDI_type = str(getattr(row, 'DDI_type'))
        sentence = getattr(row, 'Sentence')
        score = getattr(row, 'Score')
        std = getattr(row, 'STD')
        Confidence_DDI = 0
        left_drug_side_effect = ''
        right_drug_side_effect = ''
        if left_drug in drug_side_effect_info:
            left_drug_side_effect = drug_side_effect_info[left_drug]
        if right_drug in drug_side_effect_info:
            right_drug_side_effect = drug_side_effect_info[right_drug]

        if score-std/2 > model_threshold:
            Confidence_DDI = 1
            
        left_corresponding_drugs = left_ddi_info[DDI_type]
        right_corresponding_drugs = right_ddi_info[DDI_type]
        
        left_drug_similarity_df = similarity_df.loc[left_drug][left_corresponding_drugs]
        left_selected_drugs = list(left_drug_similarity_df[left_drug_similarity_df>=structure_threshold].index)
        
        right_drug_similarity_df = similarity_df.loc[right_drug][right_corresponding_drugs]
        right_selected_drugs = list(right_drug_similarity_df[right_drug_similarity_df>=structure_threshold].index)

        left_drug_annotation_string = ';'.join(left_selected_drugs)
        right_drug_annotation_string = ';'.join(right_selected_drugs)
        fp.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(prescription, drug_pair, DDI_type, sentence, 
        score, std, Confidence_DDI, left_drug_side_effect, right_drug_side_effect,
         left_drug_annotation_string, right_drug_annotation_string, left_drug, right_drug))
    
    e = time.time()

    fp.close()        
    return

def summarize_prediction_outcome(result_file, output_file, information_file):    
    sentence_interaction_info = read_information_file(information_file)
    
    with open(result_file, 'r') as fp:
        fp.readline()
        out_fp = open(output_file, 'w')
        out_fp.write('%s\t%s\t%s\t%s\t%s\t%s\n'%('Prescription', 'Drug_pair', 'DDI_type', 'Sentence', 'Score', 'STD'))
        for line in fp:
            sptlist = line.strip().split('\t')
            drug_pair_info = sptlist[0].strip()
            drug_pair_list = drug_pair_info.split('_')
            prescription = drug_pair_list[0]
            drug1 = drug_pair_list[1]
            drug2 = drug_pair_list[2]
            drug_pair = '%s_%s'%(drug1, drug2)
            DDI_class = sptlist[1].strip()
            predicted_score = sptlist[2].strip()
            predicted_std = sptlist[3].strip()
        
            template_sentence = sentence_interaction_info[DDI_class]
            prediction_outcome = template_sentence.replace('#Drug1', drug1)
            prediction_outcome = prediction_outcome.replace('#Drug2', drug2)
            out_fp.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(prescription, drug_pair, DDI_class, prediction_outcome, predicted_score, predicted_std))
        out_fp.close()

def processing_network(df, type_df):    
    PK_type_list = list(type_df['type'])    
    type_to_action = {}
    type_to_perpetrator = {}
    for row in type_df.itertuples():
        type_num = getattr(row, 'type')
        action = getattr(row, 'action')
        perpet = getattr(row, 'perpetrator')
        type_to_action[type_num] = action
        type_to_perpetrator[type_num] = perpet    
    PK_df= df[df['Interaction_type'].isin(PK_type_list)]
    action_list = []
    perpet_list = []
    victim_list = []    
    for row in PK_df.itertuples():
        type_num = getattr(row, 'Interaction_type')
        drug1 = getattr(row, 'drug1')
        drug2 = getattr(row, 'drug2')
        action = type_to_action[type_num]
        perpet = type_to_perpetrator[type_num]
        action_list.append(action)  
        if perpet == '#Drug2':
            perpet_list.append(drug2)
            victim_list.append(drug1)
        else:
            perpet_list.append(drug1)
            victim_list.append(drug2)            
    PK_df['action'] = action_list        
    PK_df['perpetrator'] = perpet_list
    PK_df['victim'] = victim_list
    return PK_df

def _get_unidirectional_pred(tmp):
    direction = None
    max_key = None
    standard = -float('inf')
    for key, val in tmp.items():
        if val[1] > standard:
            standard = val[1]
            max_key = key
            direction = val[0]
            
        elif val[1] == standard and direction != val[0]:
            max_key = None

    if max_key == None:
        return {}
    else:
        return {key: val for key, val in tmp.items() if val[0]==direction}

def find_conflicts(df):

    reported_double = {}
    reported_case = []
    for row in df.itertuples():
        drug_pair_in_order = (row[3], row[4])
        if drug_pair_in_order not in reported_double:
            reported_double[drug_pair_in_order] = (row[0], row[-2])
        elif row[-2] != reported_double[drug_pair_in_order][1]:
            reported_case += [row[0], reported_double[drug_pair_in_order][0]]
            
    df = df.loc[reported_case]
    
    severity_score_dict = {'Major': 5, 'Moderate':4, 'Minor':3, 'Not severe':2, 'Unknown':1}
    conflicting_pairs = {}
    for row in df.itertuples():
        perpet = getattr(row, 'perpetrator')
        victim = getattr(row, 'victim')
        pair = (perpet,victim)
        if pair not in conflicting_pairs:
            conflicting_pairs[pair] = {}
            conflicting_pairs[pair][row[0]] = getattr(row, 'action'), severity_score_dict[getattr(row, 'Severity')]
        else:
            conflicting_pairs[pair][row[0]] = getattr(row, 'action'), severity_score_dict[getattr(row, 'Severity')]
    idx_list = []
    for drug_pair, conflicted_pred in conflicting_pairs.items():
        filtered_result = _get_unidirectional_pred(conflicted_pred)
        if len(filtered_result) > 0:
            for k, v in filtered_result.items():
                idx_list.append(k)
    final = list(set(reported_case)-set(idx_list))
    return final

def filter_final_result(annotated_result_file, conflicting_type_file, output_file):    
    df = pd.read_csv(annotated_result_file, sep = '\t')
# Filter 1: DDI prediction confidence == 1 
    df1 = df[df['Confidence_DDI']==1]
# Filter 2: at least one structural analog with the same DDI type in gold standard dataset
    df1 = df1[(df1['Similar approved drugs (left)'].isna() == False)|(df1['Similar approved drugs (right)'].isna() == False)]

# Filter 3: severity -- 'not severe' filter out  --> eliminated filter
    # df1 = df1[df1['Final severity'] != 'Not severe']
    df1 = df1[['drug1', 'drug2', 'Interaction_type', 'Sentence', 'Final severity', 'Side effects (left)', 'Side effects (right)']]
    df1.rename(columns = {'Final severity':'Severity'}, inplace = True)

# Filter 4: drop duplicates for the model2 types with no direction (#drug1-#drug2)
    false_types = [117, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 143, 145, 146, 147, 148, 151, 153, 156,
     157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 179, 186, 187, 188, 189, 190, 191, 192, 193,
      194, 195, 196, 197, 198, 199, 200]

    pairset_list = []
    for row in df1.itertuples():
        drug1 = getattr(row, 'drug1')
        drug2 = getattr(row, 'drug2')
        pairset = frozenset([drug1, drug2])
        pairset_list.append(pairset)
    df1['pairset'] = pairset_list
    df_sub1 = df1[df1['Interaction_type'].isin(false_types)]
    df_sub1.drop_duplicates(subset = ['Interaction_type', 'pairset'], inplace = True)
    df_sub2 = df1[~df1['Interaction_type'].isin(false_types)]
    df1 = pd.concat([df_sub1, df_sub2], ignore_index = True)

# Filter 5: filter out drug pairs that have conflicting DDI types
    conflicting_types = pd.read_csv(conflicting_type_file, sep = '\t')
    conc_type_df = conflicting_types[conflicting_types['category']=='concentration']
    metab_type_df = conflicting_types[conflicting_types['category']=='metabolism']
    Conc = processing_network(df1, conc_type_df)
    metabolism = processing_network(df1, metab_type_df)
    Conc_reported = find_conflicts(Conc)
    met_reported = find_conflicts(metabolism)
    conflicts_total = met_reported + Conc_reported
    no_conflicts = list(set(df1.index)-set(conflicts_total))
    df_final = df1.loc[no_conflicts]
    df_final = df_final[['drug1', 'drug2', 'Interaction_type', 'Sentence', 'Severity', 'Side effects (left)', 'Side effects (right)']]   
    df_final.to_csv(output_file, sep = '\t', index = False)

def concatenate_results(model1_result_file, model2_result_file, output_file):
    df_model1 = pd.read_csv(model1_result_file, sep = '\t')
    df_model2 = pd.read_csv(model2_result_file, sep = '\t')
    df_concat = pd.concat([df_model1, df_model2], ignore_index = True)
    df_concat = df_concat.drop_duplicates()
    df_concat.to_csv(output_file, sep = '\t', index = False)
    
