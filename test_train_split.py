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
	
	tr_atoms = atoms[atoms['molecule_name'].isin(tr_names)]
	tr_atoms.reset_index(inplace=True, drop=True)

	te_atoms = atoms[atoms['molecule_name'].isin(te_names)]
	te_atoms.reset_index(inplace=True, drop=True)

	print('pairs..')
	
	tr_pairs = pairs[pairs['molecule_name'].isin(tr_names)]
	tr_pairs.reset_index(inplace=True, drop=True)
	
	te_pairs = pairs[pairs['molecule_name'].isin(te_names)]
	te_pairs.reset_index(inplace=True, drop=True)

	return tr_atoms, tr_pairs, te_atoms, te_pairs
