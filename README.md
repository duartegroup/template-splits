# template-splits

Code to complement our paper: An exploration of dataset importance in single-step retrosynthesis prediction. 

This repository contains the code for cleaning and splitting USPTO and Pistachio datasets based on templates, as well as the configuration files used to train models used in the paper.

## 1. Installation/Setup

1. Clone the repository:
    `git clone git@github.com:saratanov/template-splits.git`
2. Install the required Python environment using the provided yaml file:  
    `conda env create -f envs/environment.yml`
3. Activate environment:
    `conda activate template-splits`.
4. Add the directory of this repo to your Python path:    
    `export PYTHONPATH=${PYTHONPATH}:$(pwd)` (from inside the directory this README lives in).

## 2. Preprocessing

Note that we take advantage of high memory availability to load and process each dataset in full. Use the --uspto-only flag on each script if not using the Pistachio dataset.

### 2.1 Initialisation

As Pistachio is a proprietary database owned by NextMove Software, we do not provide the data here. Instead it must be manually imported into `preprocessing/data/raw` as `pistachio.smi` before running the following script. 

`1_initialisation.py`:  Initialises USPTO and Pistachio datasets for further preprocessing. The Pistachio dataset is filtered to remove US patents and remapped using RXNMapper (GPU recommended).

### 2.2 Cleaning

`2_cleaning.py`: Cleans and filters atom-mapped reaction strings for the initialised USPTO and Pistachio datasets according to the steps detailed in S1.

### 2.3 Template filtering

Templates must be extracted from the intermediate datasets `data/raw/2_uspto.csv` and `data/raw/2_pistachio.csv` prior to running the following script. We do this using the LocalTemplate module of LocalRetro: https://github.com/kaist-amsg/LocalRetro. This can also be done with alternative template extractors, such as RDChiral: https://github.com/connorcoley/rdchiral. Template SMARTS must be stored in a column named 'template'.

`3_template_filtering.py`: Filters USPTO and Pistachio by the number of template occurrences, removing all templates with less than 5 or 20 reactions respectively.

### 2.4 Template splitting

`4_template_splitting.py`: Splits the USPTO dataset via the narrow and broad template splitting strategies, as detailed in the paper. Splits the Pistachio dataset into two 10k reaction test sets, where the ID and OOD test sets contains templates present or not present in USPTO respectively.

#### 3. Training

All configuration files used to train the models in the paper are provided. Models were trained according to their respective packages:
    LocalRetro: https://github.com/kaist-amsg/LocalRetro (since removed by the authors)
    MEGAN: https://github.com/molecule-one/megan 
    RootAligned: https://github.com/otori-bird/retrosynthesis