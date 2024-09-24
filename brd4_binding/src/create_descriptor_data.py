"""
Create descriptors file from data with BuildingBlock SMILES.

Read parquet file with pandas and creates and saves a
parquet file with the same ids from original file and all
chemical descriptors from RDKit as columns with the values
from each BuildingBlock in a 3-tuple.

"""

import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors


original_df = pd.read_parquet("../data/brd4_5p.parquet").set_index('id')

# Remove protecting group from buildingblock1
prot_group = 'CC1c2ccccc2-c2ccccc21'
for index in original_df.index:
    old_str = original_df.loc[index, 'buildingblock1_smiles']
    new_val = old_str.replace(prot_group, '')
    original_df.at[index, 'buildingblock1_smiles'] = new_val

# Get descriptors from Chem.Descriptors
descriptors_list = [x[0] for x
                    in Descriptors._descList
                    ]
total_descriptors = len(descriptors_list)


# Create dataframes with descriptors from buildingblocks
def make_descriptor_df(buildingblock):
    bb_smiles = {y: x for x, y in
                 enumerate(original_df["buildingblock{}_smiles"
                                       .format(buildingblock)].unique())}
    descriptors = [Descriptors.CalcMolDescriptors(Chem.MolFromSmiles(smile))
                   for smile in list(bb_smiles.keys())]
    descriptors_df = pd.DataFrame(descriptors)
    index = "bb{}_smiles".format(buildingblock)
    descriptors_df[index] = list(bb_smiles.keys())
    descriptors_df = descriptors_df.set_index(index)

    return descriptors_df


descriptors_bb1_df = make_descriptor_df(1)
descriptors_bb2_df = make_descriptor_df(2)
descriptors_bb3_df = make_descriptor_df(3)
