# **Data Processing**

## _This repository contains scripts that help in the preprocessing of sdf files and molecular dataframes before making an IMPRESSION prediction or analysing the predictions afterwards_

### Requirements:
#### - RDKit
#### - Pandas
#### - Numpy


## *Conformer Generation*

#### This script will sample molecules from a dataset and generate a certain number of conformers for each molecule (via RDKit) depending on the number of rotational bonds in the molecule and then perform a redundant elimination to remove any conformers that are too similar to eachother in geometry and energy.

## *Checks*

#### This script contains functions for running checks on molecules before submitting them to IMPRESSION. The checks involve ensuring only molecules with the compatible nuclei* are included and 
