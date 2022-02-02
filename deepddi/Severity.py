import os
import glob
import pickle
import copy
import argparse
import numpy as np
import time

import pandas as pd
from keras.models import model_from_json 
from sklearn.preprocessing import MultiLabelBinarizer

def predict_severity(output_dir, output_file, df, trained_model, trained_weight, threshold):  
    with open('./data/multilabelbinarizer_severity.pkl', 'rb') as fid:
        lb = pickle.load(fid)
    
    #df = pd.read_csv(ddi_input_file, index_col=0)    
    ddi_pairs = list(df.index)
    X = df.values    
    
    json_file = open(trained_model, "r")
    loaded_model_json = json_file.read() 
    json_file.close()
    
    model = model_from_json(loaded_model_json)
    model.load_weights(trained_weight)
    
    mc_predictions = []
    iter_num = 10
    for i in range(iter_num):
        y_predicted = model.predict(X)
        mc_predictions.append(y_predicted)

    arr = np.asarray(mc_predictions)
    
    y_predicted_mean = np.mean(arr, axis=0)
    y_predicted_std = np.std(arr, axis=0)
    
    original_predicted_ddi = copy.deepcopy(y_predicted_mean)
    original_predicted_ddi_std = copy.deepcopy(y_predicted_std)
    
    y_predicted_mean[y_predicted_mean >= threshold] = 1
    y_predicted_mean[y_predicted_mean < threshold] = 0
    
    y_predicted_inverse = lb.inverse_transform(y_predicted_mean)   
    
    fp = open(output_file, 'w')
    fp.write('Drug pair\tPredicted class\tScore\tSTD\n')
    for i in range(len(ddi_pairs)):
        predicted_ddi_score = original_predicted_ddi[i]

        predicted_ddi_std = original_predicted_ddi_std[i]
        predicted_ddi = y_predicted_inverse[i]
        each_ddi = ddi_pairs[i]           
        for each_predicted_ddi in predicted_ddi:
            fp.write('%s\t%s\t%s\t%s\n'%(each_ddi, each_predicted_ddi, predicted_ddi_score[each_predicted_ddi-1], predicted_ddi_std[each_predicted_ddi-1]))
    fp.close()
