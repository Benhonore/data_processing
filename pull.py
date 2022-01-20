import pandas as pd
import numpy as np
import pickle
import os
import glob
import shutil
from rdkit import Chem
from rdkit.Chem import AllChem


# This function is used to identify molecules containing only a certain set of nuclei within it. It takes an input file in nmredata format. 
# It returns false if the molecule containes nuclei outside of the specified list and true if only nuclei from the list are present.
# It is useful to loop through a set of molecules using this function.

def pull_nuclei_mols(file, atoms = []):
  
  nuclei = []
  if 'H' in atoms:
    nuclei.append(1)
  if 'C' in atoms:
    nuclei.append(6)
  if 'N' in atoms:
    nuclei.append(7)
  if 'O' in atoms:
    nuclei.append(8)
  if 'F' in atoms:
    nuclei.append(9)
  if 'Si' in atoms:
    nuclei.append(14)
  if 'P' in atoms:
    nuclei.append(15)
  if 'S' in atoms:
    nuclei.append(16)
  if 'Cl' in atoms:
    nuclei.append(17)
  if 'Br' in atoms:
    nuclei.append(35)
    
  atomic_nums = []
  
  with open(file, 'r') as f:
    
    for line in f:
      
      if len(line.split(',')) == 4:
          
          atommic_num = int(line.split(',')[2])
          atomic_nums.append(atomic_num)
      
  if not all([x in nuclei for x in atomic_nums]):
    
    return 'False'
  
  else:
    
    return 'True'
             
  if len(atomic_nums) == 0:
    print('Nuclei are not being read properly. What file type are you using? This method reads nmredata sdf files.')        
        
        
        
  
