import pandas as pd
import numpy as np
from tqdm import tqdm

### This function reads in normal dataframes made from 3D molecular structures and turns them into 2D by averaging 
### 2D-equivalent protons, changing distance to bond order and zeroing out coupling/nmr_type/path_len


def make_2d(atom_df, pair_df, name):

	print('updating atoms..')

	df=pd.read_pickle(atom_df)
	df['averaged_shift'] = 0

	for molecule in tqdm(np.unique(df['molecule_name'])):

		mol_df = df[df['molecule_name'] == molecule]
		d = {}

		for atom in range(len(mol_df)):

			if mol_df.iloc[atom]['typestr'] == 'C' or mol_df.iloc[atom]['typestr'] == 'N':
				protons = {}

				for i, con in enumerate(mol_df.iloc[atom]['conn']):
					if con != 0 and mol_df.iloc[i]['typestr'] == 'H':
						protons[i] = float(mol_df.iloc[i]['shift'])

				if len(protons) !=0:
					d[atom]=protons

		ave = {}

		for j in d:
			for h in d[j]:
				ave[h] = np.mean(list(d[j].values()))

		for proton in ave:
			i = mol_df.index[mol_df['atom_index']==proton].to_list()[0]
			mol_df.at[i, 'averaged_shift'] = ave[proton]

		df[df['molecule_name']==molecule] = mol_df


	for i in range(len(df)):
		if df.iloc[i]['averaged_shift'] ==0:
			df.at[i, 'averaged_shift'] = df.iloc[i]['shift']

	df.rename(columns={'shift':'3d_shift', 'averaged_shift':'shift'}, inplace=True)
	
	df.to_pickle('2d_' + name + '_atoms.pkl')

	print('updating pairs..')

	bo = []
	
	for atom in tqdm(range(len(df))):
		for i in df.iloc[atom]['conn']:
			bo.append(i)

	dc = pd.read_pickle(pair_df)
	
	indices = ['nmr_type', 'coupling', 'path_len']
	for i in indices:
		dc[i] = 0

	dc['3d_distance'] = dc['dist']
	
	assert len(dc) == len(bo), 'problem! atoms and pairs dont seem to line up..'

	dc['dist'] = bo

	dc.to_pickle('2d_' + name + '_pairs.pkl')

