import argparse
import os
import shutil
import logging
import time

from deepddi import DeepDDI
from deepddi import Severity
from deepddi import preprocessing
from deepddi import result_processing

if __name__ == '__main__':
    start = time.time()
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-o', '--output_dir', required=True, help="Output directory")
    parser.add_argument('-i', '--input_file', required=True, help="Input file")
    parser.add_argument('-p', '--PCA_profile_file', help="PCA profile file")
    # parser.add_argument('-m', '--ddi_trained_model', required=True, help = "'drugbnak' or 'manual'")
    
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
    
    options = parser.parse_args()
    raw_input_file = options.input_file
    output_dir = options.output_dir
    PCA_profile_file = options.PCA_profile_file
    
    drug_dir = './data/DrugBank5.0_Approved_drugs/'
    pca_model = './data/PCA_tanimoto_model_50.pkl'
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    drug_information_file = './data/Approved_drug_Information.txt'
    drug_enzyme_information_file = './data/Approved_drug_enzyme_Information.txt'   
    side_effect_information_file = './data/Drug_Side_Effect.txt'
    
    drug_list = []
    with open('./data/DrugList.txt', 'r') as fp:
        for line in fp:
            drug_list.append(line.strip())

    input_file = '%s/DDI_input.txt' % (output_dir)
    preprocessing.parse_DDI_input_file(raw_input_file, input_file)
    
    known_drug_similarity_file = './data/drug_similarity.csv'
    similarity_profile = known_drug_similarity_file


    if PCA_profile_file == None:
        similarity_profile = '%s/similarity_profile.csv' % output_dir
        pca_similarity_profile = '%s/PCA_transformed_similarity_profile.csv' % output_dir
        pca_profile_file = '%s/PCA_transformed_similarity_profile.csv' % output_dir
        print ('Calculate structure similarity profile')
        preprocessing.calculate_structure_similarity(drug_dir, input_file, similarity_profile, drug_list)
        preprocessing.calculate_pca(similarity_profile, pca_similarity_profile, pca_model)
        print ('Combine structural similarity profile')
        pca_df = preprocessing.generate_input_profile(input_file, pca_similarity_profile)
    else:
        pca_similarity_profile = PCA_profile_file
        print ('Generate input profile')
        pca_df = preprocessing.generate_input_profile(input_file, pca_similarity_profile)

    # threshold = 0.5
    model1_threshold = 0.4

    # model processing
    print('model running')

    model1_dir = output_dir

    ddi_trained_model = './data/models/ddi_model.json'
    ddi_trained_model_weight = './data/models/ddi_model.h5'
    DDI_sentence_information_file = './data/Type_information/Interaction_information_model1.csv'
    binarizer_file = './data/multilabelbinarizer.pkl'
    known_DDI_file = './data/DrugBank_known_ddi.txt'

    output_file = '%s/DDI_result.txt' % (model1_dir)
    ddi_output_file = '%s/Final_DDI_result.txt' % (model1_dir)   
    annotation_output_file = '%s/Final_annotated_DDI_result.txt' % (model1_dir)
    model_type = 'model1'
    DeepDDI.predict_DDI(output_file, pca_df, ddi_trained_model, ddi_trained_model_weight, model1_threshold, binarizer_file, model_type)    
    result_processing.summarize_prediction_outcome(output_file, ddi_output_file, DDI_sentence_information_file)
    result_processing.annotate_DDI_results(ddi_output_file, drug_information_file, drug_enzyme_information_file,
    similarity_profile, known_DDI_file, annotation_output_file, side_effect_information_file, model1_threshold, 0.7)    

    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))
