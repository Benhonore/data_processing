import pandas as pd
import numpy as np
from tqdm import tqdm
import random


def test_train_split(atoms, pairs, split=8):
	
	names = np.unique(atoms['molecule_name'])
	random.shuffle(names)
	
	tr_names = names[:int(float('0.'+str(split))*len(names))]
	te_names = names[int(float('0.'+str(split))*len(names)):]
	
	print('atoms..')
	atoms['filter']=0

	for i in tqdm(range(len(atoms))):
		if atoms.iloc[i]['molecule_name'] in te_names:
			atoms.at[i, 'filter'] = 1

	print('pairs..')
	pairs['filter'] = 0
	
	for i in tqdm(range(len(pairs))):
		if pairs.iloc[i]['molecule_name'] in te_names:
			pairs.at[i, 'filter'] = 1

	tr_atoms = atoms[atoms['filter'] == 0]
	tr_atoms.drop(columns=['filter'], inplace=True)
	tr_atoms.reset_index(inplace=True, drop=True)

	te_atoms = atoms[atoms['filter'] == 1]
	te_atoms.drop(columns=['filter'], inplace=True)
	te_atoms.reset_index(inplace=True, drop=True)

	tr_pairs = pairs[pairs['filter'] == 0]
	tr_pairs.drop(columns=['filter'], inplace=True)
	tr_pairs.reset_index(inplace=True, drop=True)

	te_pairs = pairs[pairs['filter'] == 1]
	te_pairs.drop(columns=['filter'], inplace=True)
	te_pairs.reset_index(inplace=True, drop=True)

	return tr_atoms, tr_pairs, te_atoms, te_pairs
