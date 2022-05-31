import re
import pandas as pd
import numpy as np
import os
import subprocess
from collections import defaultdict
# DDI

MODEL_DIR = '/Documents/GitHub/deepddi2/'
INPUT_PATH = './DDI_input.txt'
OUTPUT_DIR = './output'
OUTPUT_TXT = 'output/Final_annotated_DDI_result.txt'
SIGNIFICANCE = 0.8
DFI_INPUT_DRUGS = []
DDI_OTHER_DRUGS = []

food_comp = pd.read_csv('./database/food_compounds_lookup.csv')
# pd.DataFrame({'Name': food.orig_food_common_name\
#                           .str.replace('(', '').str.replace(')', '')\
#                           .str.split()\
#                           .str[0].unique()}
#             ).to_csv('./database/food.csv', index = None)

food_loc = food_comp.orig_food_common_name\
            .str.replace('(', '')\
            .str.replace(')', '')\
            .str.split().str[0]
comp_loc = food_comp.name.str.lower()

food2comp = defaultdict(set)
comp2food = defaultdict(set)
for f,c in zip(food_loc, comp_loc):
    food2comp[f].add(c)
    comp2food[c].add(f)
    

def regex_search(desc, pools):
    # assume desc is lowercased
    out = []
    for elem in pools:
        pattern = elem.strip().lower()
        if re.search(pattern, desc):
            out.append(elem)
    return out

def ingest_input(input_json, interaction_type, input_fp = INPUT_PATH,
                 compounds_path = './database/drug_info_combined.csv',
                 food_path = './database/food.csv'):
    assert interaction_type.lower() in ['ddi', 'dfi'], 'API not supported'
    first_line = []
    second_line = []
    # handle ddi
    if interaction_type.lower() == 'ddi':
        drug_pools = pd.read_csv(compounds_path).Name.str.lower()
        
        cur_desc = input_json['current_drug']['drug_desc'].lower()
        drug_title = input_json['current_drug']['drug_title']
        
        drug_search = regex_search(cur_desc.lower(), drug_pools)
        assert drug_search, ('Drug: %s not Found' % drug_title)
            
        first_line += [drug_title+'|'+i for i in drug_search]
        for drug in input_json['other_drug']:
            DDI_OTHER_DRUGS.append(drug['drug_title'].lower())
            other_desc = drug['drug_desc'].lower()
            other_search = regex_search(other_desc, drug_pools)
            
            if not other_search:
                print('Drug: %s not Found' % drug['drug_title'])
                continue
            second_line += [drug['drug_title']+'|'+ i for i in other_search]
    # handle dfi
    else:
        drug_pools = pd.read_csv(compounds_path).Name.str.lower()
        food_pools = pd.read_csv(food_path).Name
        for drug in input_json['drug_list']:
            DFI_INPUT_DRUGS.append(drug['drug_title'].lower())
            drug_search = regex_search(drug['drug_desc'].lower(), drug_pools)
#             print(drug['drug_desc'], drug_search)
            if not drug_search:
                print('Drug: %s not Found' % drug['drug_title'])
                continue
            first_line += [drug['drug_title'] + '|' + i for i in drug_search]
        
        for food in input_json['food_list']:
            food_search = food2comp[food] 
            if not food_search:
                print('Food: %s not Found' % food)
                continue
            second_line += [food + '|' + i for i in food_search]
            
    # TODO handle not-found case           
    if os.path.exists(input_fp):
        os.remove(input_fp)
    
    with open(input_fp, 'w') as fw:
        fr = '\t'.join(first_line) + '\n'
        sr = '\t'.join(second_line) + '\n'
        fw.write(fr)
        fw.write(sr)


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


def collect_output(thres = SIGNIFICANCE, out_txt = OUTPUT_TXT):
    res = pd.read_csv(out_txt,
                      sep='\t', 
                      header=0)[['drug1', 'drug2',
                                          'DDI_prob', 'DDI_prob_std',
                                          'Confidence_DDI', 'Sentence',
                                          'Side effects (left)',
                                          'Side effects (right)']]
    
    temp = res.loc[(res.Confidence_DDI == 1) &\
                   (res.DDI_prob >= thres)]
    
    if DDI_OTHER_DRUGS:
        temp = temp.loc[(res.drug1.str.contains(r'|'.join(DDI_OTHER_DRUGS))) |\
                        (res.drug2.str.contains(r'|'.join(DDI_OTHER_DRUGS)))]
        
    if DFI_INPUT_DRUGS:
        temp = temp.loc[(res.drug1.str.contains(r'|'.join(DFI_INPUT_DRUGS))) |\
                        (res.drug2.str.contains(r'|'.join(DFI_INPUT_DRUGS)))]
    
    out = {}
    out['interactions'] = []
    for line in temp.iterrows():
        inner_out = {}
        row = line[1].values
        other_drug = row[1] ###
        interaction_desc = row[5]
        side_effect = {}
        side_effect['probability'] = row[2]
        side_effect['side_effect_id'] = row[6]
        inner_out['other_drug'] = other_drug
        inner_out['interaction_desc'] = interaction_desc
        inner_out['side_effect'] = side_effect
        out['interactions'].append(inner_out)
    return out

# Example of calling DFI API
dfi_sample_input = {'drug_list': [{'drug_title': 'Drug C', 'drug_desc': '? Aspirin '}, {'drug_title': 'Drug D', 'drug_desc': '? Vitamin C '}],
 'food_list': ['lemon', 'orange']}

# Example of calling DDI API
ddi_sample_input =  {'current_drug': {'drug_title': 'Good Drug', 'drug_desc': 'Vitamin C'},
 'other_drug': [{'drug_title': 'Drug A', 'drug_desc': ' cool  Vitamin A '},
                {'drug_title': 'Drug B', 'drug_desc': ' very good Aspirin Acetaminophen'}]}

# ALL you need to call is func 'run(input_json,type)'
print(run(ddi_sample_input,'DDI'))
# print(run(ddi_sample_input,'DDI'))

    