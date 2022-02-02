# DeepDDI-X #
#Project
This project is to develop a framework that systematically predicts drug-drug interactions (DDIs).

#Features
- DeepDDI-X predicts DDIs using names of drug-drug pairs and their structural information (SMILES) as inputs.
- DeepDDI-X predicts 201 DDI types in human-readable sentences and 4 levels of clinical severity as output.

###Procedure
**Note**: This source code was developed in Linux, and has been tested in Ubuntu 16.04

1. Clone the repository

        git clone https://github.com/kaistsystemsbiology/DeepDDIX/

2. Create and activate a conda environment

        conda env create -f environment.yml
        conda activate deepddi

##Input files for the DeepDDI##
Following working input files can be found in: `./data/`. These files were used for the data presented in the manuscript.

#Example#
```
python run_DeepDDI.py -i ./DDI_input.txt -o ./output
```

