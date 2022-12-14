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

def pull_nuclei_mols(file, atoms = [], path):
        
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
          
                                atomic_num = int(line.split(',')[2])
                                atomic_nums.append(atomic_num)
                        
        Assert len(atomic_nums) != 0, 'Nuclei are not being read properly. What file type are you using? This method reads nmredata sdf files.'
      
        if not all([x in nuclei for x in atomic_nums]):
    
                print(f'Non {atoms} nuclei')
                continue
  
        else:
    
                shutil.copy(file, path)
             
       

# This function fixes a discrepancy between different nmredata file types so that RDKit can read them properly.
   
def m_end(file):
        # Read in the file   
        print(file)
        f = open(file, 'r')
        read_file = f.readlines()
        c=0
        
        # Locaste the line that M END is on
        for line in read_file:
                c+=1
        if 'END' in line:
                print(line)
                print(c-1)
                break
        print(c-1)
        print(read_file[c-1])
        
        # Amend the M END line
        read_file[c-1] = 'M  END'
        print(read_file[c-1])
        
        # Write the correct version back into the file
        g = open(file, 'w')
        g.writelines(read_file)

        
        
# For RDKit made molecules, this function pulls out any molecules that are in 2D

def find_2d_mols(file):
        with open(file, 'r') as f:
                for line in f:
                        if len(line.split()) == 2 and line.split()[0] == 'RDKit':
                                dimensionality = line.split()[1]
                                if dimensionality == '2D':
                                        print('molecule is 2D')
                                          

  
