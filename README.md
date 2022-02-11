# DeepDDI2 #

#Project
This project is to develop a framework that systematically predicts drug-drug interactions (DDIs).

#Features
- DeepDDI2 predicts DDIs using names of drug-drug pairs and their structural information (SMILES) as inputs.
- DeepDDI2 predicts 113 DDI types in human-readable sentences as output.

###Procedure
**Note**: This source code was developed in Linux, and has been tested in Ubuntu 16.04

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

