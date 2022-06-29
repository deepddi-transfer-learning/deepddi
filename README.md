# DeepDDI2 #


#Project Citation: [10.1073/pnas.1803294115](https://www.pnas.org/doi/full/10.1073/pnas.1803294115)

#Code Reference: [deepddi](https://bitbucket.org/kaistsystemsbiology/deepddi/src/master/)

#Project
This project is to develop a framework that systematically predicts drug-drug interactions (DDIs) and drug-food ineractions(DFIs). This frame work is heavily built on top of https://bitbucket.org/kaistsystemsbiology/deepddi/src/master/

#Features
- DeepDDI predicts DDIs using a list of drugs information (description + name) in json format. 
- DeepDDI predicts DFIs using single drug input (description + name) and a list of food names in json format.

1. Clone the repository

        git clone https://bitbucket.org/kaistsystemsbiology/deepddi2.git

2. Create and activate a conda environment

        conda env create -f environment.yml
        conda activate deepddi

##Input files for the DeepDDI2##
Following working input files can be found in: `./data/`. These files were used for the data presented in the manuscript.

#Example#
```
python run_DeepDDI.py -i ./DDI_input.txt -o ./output
```

Further details can be found in Supplementary Information.
