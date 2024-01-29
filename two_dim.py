import pandas as pd
import numpy as np
from tqdm import tqdm

### This function reads in normal dataframes made from 3D molecular structures and turns them into 2D by averaging 
### 2D-equivalent protons, changing distance to bond order and zeroing out coupling/nmr_type/path_len


def make_2d(atom_df, pair_df=pd.DataFrame()):

	print('updating atoms..')

	atom_df['averaged_shift'] = 0

	for molecule in tqdm(np.unique(atom_df['molecule_name'])):

		mol_df = atom_df[atom_df['molecule_name'] == molecule]
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

		atom_df[atom_df['molecule_name']==molecule] = mol_df


	for i in range(len(atom_df)):
		if atom_df.iloc[i]['averaged_shift'] ==0:
			atom_df.at[i, 'averaged_shift'] = atom_df.iloc[i]['shift']

	atom_df.rename(columns={'shift':'3d_shift', 'averaged_shift':'shift'}, inplace=True)
		

	print('updating pairs..')

	if len(pair_df)==0:
		return atom_df, pair_df

	else:	 
		bo = []	
		for atom in tqdm(range(len(atom_df))):
			bo.extend(list(atom_df.iloc[atom]['conn']))
		
		assert len(pair_df) == len(bo), 'problem! atoms and pairs dont seem to line up..'
	
		for i in ['nmr_type', 'coupling', 'path_len']:
			pair_df[i] = 0

		pair_df['3d_distance'] = pair_df['dist']
		pair_df['dist'] = bo
	
		return atom_df, pair_df
