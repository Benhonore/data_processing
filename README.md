# **Data Processing**

## _This repository contains scripts that help in the preprocessing of sdf files and molecular dataframes before making an IMPRESSION prediction or analysing the predictions afterwards. It is written to compliment mol_translator._

### Requirements:

#### - mol_translator
#### - RDKit
#### - Openbabel
#### - Pandas
#### - Numpy


## *Conformer Generation*

#### This script will sample molecules from a dataset and generate a certain number of conformers for each molecule (via RDKit) depending on the number of rotational bonds in the molecule and then perform a redundant elimination to remove any conformers that are too similar to eachother in geometry and energy.

## *Checks*

#### This script contains functions for running checks on molecules before submitting them to IMPRESSION. The checks involve ensuring only molecules with the compatible nuclei* are included and that molecules are 3D or else it will embed them in 3D via RDKit.

## *Model Performance*

#### These functions 1) return a summary of the model performance at predicting carbon and proton chemical shifts and 2) in the case of conformer discrimination, will rank the closeness in prediction of the NMR carbon or proton chemical shifts for a given to the ground truth of every conformer for that molecule and return the number of conformers were ranked top.


##### * compatible nuclei: H, C, N, O, F, Si, P, S, Cl, Br


