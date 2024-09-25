"""
Module to create pandas.DataFrame of descritors for each buildingblock.

Module to create pandas.DataFrame with all buildingblocks as index and
RDKit's 2D Descriptors as columns.
"""

import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors


def _get_smiles_list(data, buildingblock) -> list:

    col = "buildingblock{}_smiles".format(buildingblock)
    return data[col].unique().tolist()


def get_descriptors() -> list:
    descriptors_list = [x[0] for x
                        in Descriptors._descList
                        ]
    return descriptors_list


# Create dataframes with descriptors from buildingblocks
def make_descriptor_df(original_data,
                       buildingblock,
                       ) -> pd.DataFrame:

    print("Parsing unique buildingblock smiles...")
    smiles_list = _get_smiles_list(original_data, buildingblock)

    print("Calculating descriptors"
          " for buildingblock{} ({})..."
          .format(buildingblock, len(smiles_list))
          )
    descriptors = [Descriptors.CalcMolDescriptors(Chem.MolFromSmiles(smile))
                   for smile in smiles_list]

    print("Creating DataFrame...")
    descriptors_df = pd.DataFrame(descriptors)
    index = "bb{}_smiles".format(buildingblock)
    descriptors_df[index] = smiles_list
    descriptors_df = descriptors_df.set_index(index)
    print("Done.")

    return descriptors_df
