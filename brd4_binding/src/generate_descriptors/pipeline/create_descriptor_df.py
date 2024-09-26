"""
Module to create pandas.DataFrame of descritors for each buildingblock.

Module to create pandas.DataFrame with all buildingblocks as index and
RDKit's 2D Descriptors as columns.
"""

import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors


def _get_smiles_list(data) -> list:

    unique_buildingblocks = []

    for i in [1, 2, 3]:
        col = "buildingblock{}_smiles".format(i)
        bb_list = list(set(data.loc[:, col].to_list()))
        unique_buildingblocks += bb_list

    return unique_buildingblocks


def get_descriptors() -> list:
    descriptors_list = [x[0] for x
                        in Descriptors._descList
                        ]
    return descriptors_list


# Create dataframes with descriptors from buildingblocks
def make_descriptor_df(original_data,
                       ) -> pd.DataFrame:
    unique_bbs = _get_smiles_list(original_data)
    temp_dict = {}

    for bb in unique_bbs:
        temp_dict[bb] = Descriptors.CalcMolDescriptors(Chem.MolFromSmiles(bb))

    descriptor_df = pd.DataFrame.from_dict(temp_dict, orient="index")

    return descriptor_df
