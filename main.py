import re
import pandas as pd
import numpy as np
import os
import subprocess
# DDI

MODEL_DIR = '/Documents/GitHub/deepddi2/'
INPUT_PATH = './DDI_input.txt'
OUTPUT_DIR = './output'
OUTPUT_TXT = 'output/Final_annotated_DDI_result.txt'
SIGNIFICANCE = 0.9
DFI_INPUT_DRUGS = []

def regex_search(desc, pools):
    # assume desc is lowercased
    out = []
    for elem in pools:
        pattern = elem.strip().lower()
        if re.search(pattern, desc):
            out.append(elem)
    return out

def ingest_input(input_json, interaction_type, input_fp = INPUT_PATH,
                 compounds_path = './database/Drug_info_combined.csv',
                 food_path = './database/food.csv'):
    assert interaction_type.lower() in ['ddi', 'dfi'], 'API not supported'
    out = []
    
    # handle ddi
    if interaction_type.lower() == 'ddi':
        drug_pools = pd.read_csv(compounds_path).Name
        
        cur_desc = input_json['current_drug']['drug_desc'].lower()
        drug_name = input_json['current_drug']['drug_title']
        out += [drug_name+'|'+i for i in regex_search(cur_desc, drug_pools)]
        for drug in input_json['other_drug']:
            each_name = drug['drug_title']
            other_desc = drug['drug_desc'].lower()
            found=[each_name+'|'+ i for i in regex_search(other_desc, drug_pools)]
            out += found
    # handle dfi
    else:
        drug_pools = pd.read_csv(compounds_path).Name
        food_pools = pd.read_csv(food_path).Name
        cur_desc = input_json['current_drug']['drug_desc'].lower()
        drug_name = input_json['current_drug']['drug_title']
        out += [drug_name+'|'+i for i in regex_search(cur_desc, drug_pools)][:1]
        # handle drug not found
        assert len(regex_search(cur_desc, drug_pools))>0 , 'We can not found drug: '+drug_name
        for food in input_json['food']:
            out+= regex_search(food.lower(), food_pools)      
        
    if os.path.exists(input_fp):
        os.remove(input_fp)

    with open(input_fp, 'w') as fw:
        to_write = '\n'.join(out) + '\n'
        fw.write(to_write)


def run(input_json,interaction_type,thres=SIGNIFICANCE):
    # INPUT: 
    #   input_json: the json file of input info
    #   interactioin_type: 'DFI' or 'DDI'

    # execute & make sure it runs linearly
    ingest_input(input_json, interaction_type)
    cmd = ('python run_DeepDDI.py -i %s -o %s -t %s'%('DDI_input.txt', 'output',str(interaction_type))).split()
    try:
        subprocess.run(cmd)
        return collect_output(thres)
    except AssertionError:
        return None

def format_output():
    # TODO: handle DFI comp2food 
    # TODO: hanlde interaction id2desc
    pass

def collect_output(thres = SIGNIFICANCE, out_txt = OUTPUT_TXT):
    thres = 0.95
    res = pd.read_csv(out_txt,
                      sep='\t', 
                      header=0)[['drug1', 'drug2',
                                          'DDI_prob', 'DDI_prob_std',
                                          'Confidence_DDI', 'Sentence',
                                          'Interaction_type']]
    
    temp = res.loc[(res.Confidence_DDI == 1) &\
                   (res.DDI_prob >= thres)]
    
    if DFI_INPUT_DRUGS:
        temp = temp.loc[(res.drug1.str.contains('|'.join(DFI_INPUT_DRUGS))) |\
                        (res.drug2.str.contains('|'.join(DFI_INPUT_DRUGS)))]
    out = {}
    out['interations'] = []
    for line in temp.iterrows():
        inner_out = {}
        row = line[1].values
    #     print(row)
        drug_pair = [row[0], row[1]]
        interaction_desc = row[5]
        side_effect = {}
        side_effect['probability'] = row[2]
        side_effect['side_effect_id'] = row[6]
        inner_out['drug_pair'] = drug_pair
        inner_out['interaction_desc'] = interaction_desc
        inner_out['side_effect'] = side_effect
        out['interations'].append(inner_out)
    return out

# Example of calling DFI API
dfi_sample_input = {'current_drug': {'drug_title': 'Bad Drug', 'drug_desc': '?xxxxxx xxx Acetaminophen xxx xxx'},
 'food': ['lemon', 'orange']}

# Example of calling DDI API
ddi_sample_input =  {'current_drug': {'drug_title': 'Good Drug', 'drug_desc': '! Biotin !'},
 'other_drug': [{'drug_title': 'So-so Drug', 'drug_desc': ' cool  Vitamin C and vitamin A'},
                {'drug_title': 'Pain-Killer', 'drug_desc': ' very good Aspirin Acetaminophen'}]}

# ALL you need to call is func 'run(input_json,type)'
print(run(dfi_sample_input,'DFI'))
# print(run(ddi_sample_input,'DDI'))

    