import pandas as pd
import numpy as np

import glob
import os
import pickle
import sys
import shutil

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Lipinski

sys.path.append('INSERT MOL TRANSLATOR PATH')

from mol_translator.aemol import Aemol
from mol_translator.conformer.aeconf import Aeconf

def get_num_conformers(num_rot_bonds):
	if 3**(num_rot_bonds) < 50:
		return 3**(num_rot_bonds)
	else:
		return 50

files=glob.glob('INSERT FILE INPUT PATH')
path='INSERT FAILED PATH'

ids=[]

for file in files:
	try:
		suppl = Chem.SDMolSupplier(file, removeHs=False)
		id=file.split('/')[-1].split('.')[0]
		for mol in suppl:
			num_rot_bonds = Lipinski.NumRotatableBonds(mol)
			num_confs = get_num_conformers(num_rot_bonds)
			cids=AllChem.EmbedMultipleConfs(mol, int(num_confs))
			res=AllChem.MMFFOptimizeMoleculeConfs(mol)
			for cid in range(mol.GetNumConformers()):
				w=Chem.SDWriter(f'CREATE CONFORMERS DIR/{id}_conf_{cid}.sdf')
				w.write(mol, cid)

		with open(f'CREATE CONFORMERS DIR/conformer_energies/{id}_conformer_energies.pkl', 'wb') as fh:
			pickle.dump(res, fh)
		ids.append(id)
	except:
		shutil.move(file, path)

conf_files=glob.glob('CREATE CONFORMERS DIR/*sdf')

for id in ids:
	print(id)
	res=pd.read_pickle(f'CREATE CONFORMERS DIR/{id}_conformer_energies.pkl')
	amols=[]
	for file in conf_files:
		if file.split('/')[-1].split('_conf')[0]==id:
			cid=file.split('/')[-1].split('conf_')[-1].split('.')[0]
			amol=Aemol(f'{id}_{cid}')
			amol.from_file_ob(file, 'sdf')
			amol.from_ob(amol.obmol)
			amol.prop_from_file(file, 'nmr', 'nmredata')
			amol.mol_properties['energy']=res[int(cid)][1]
			print(amol.mol_properties)
			amols.append(amol)

	confs= Aeconf(amols)
	confs.calc_pops()
	confs.boltzmann_average()
	confs.get_dist_array()
	confs.redundant_elimination(0.35)
	
	elm_ids=[]
	
	with open('conf_eliminated_molecules.txt', 'rb') as fh:
		for line in fh.readlines():
			elm_ids.append(int(line))
	
	if len(elm_ids) != 0:
		for elm_id in elm_ids:
			os.remove(f'CREATE CONFORMERS DIR/{id}_conf_{elm_id}.sdf')

	print(f'{len(elm_ids)} molecules eliminated')
	  	
