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
from rdkit.rdBase import BlockLogs


def read_data(path="../data/brd4_5p.parquet",
              index='id',
              ) -> pd.DataFrame:
    """
    Read data from parquet file into a pandas.DataFrame.

    Read data from parquet file located in path into a
    pandas.DataFrame, sets the index to 'id' column

    Args:
        path (str): path to the parquet file
        index (str): column of the dataframe to use as index.

    Returns:
        pandas.DataFrame

    """
    return pd.read_parquet(path).set_index(index)


# Remove protecting group from buildingblock1
def clean_bb1(original_data: pd.DataFrame) -> pd.DataFrame:
    prot_group = 'CC1c2ccccc2-c2ccccc21'
    for index in original_data.index:
        old_str = original_data.loc[index, 'buildingblock1_smiles']
        new_val = old_str.replace(prot_group, '')
        original_data.at[index, 'buildingblock1_smiles'] = new_val
    return original_data


# Get descriptors from Chem.Descriptors
def get_descriptors() -> int:
    descriptors_list = [x[0] for x
                        in Descriptors._descList
                        ]
    return len(descriptors_list)


def get_smiles_list(data, buildingblock) -> list:

    col = "buildingblock{}_smiles".format(buildingblock)
    return data[col].unique().tolist()


# Create dataframes with descriptors from buildingblocks
def make_descriptor_df(original_data,
                       buildingblock,
                       ) -> pd.DataFrame:

    print("Parsing unique buildingblock smiles...")
    smiles_list = get_smiles_list(original_data, buildingblock)

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


# def get_3_tuple(id, column):
#     df_loc_bb1 = brd4_df.loc[id, 'buildingblock1_smiles']
#     bb1_value = descriptors_bb1_df.loc[df_loc_bb1, column]

#     df_loc_bb2 = brd4_df.loc[id, 'buildingblock2_smiles']
#     bb2_value = descriptors_bb2_df.loc[df_loc_bb2, column]

#     df_loc_bb3 = brd4_df.loc[id, 'buildingblock3_smiles']
#     bb3_value = descriptors_bb3_df.loc[df_loc_bb3, column]
#     return (bb1_value, bb2_value, bb3_value)


if __name__ == "__main__":
    # Silence RDKit deprecation warning
    block = BlockLogs()
    brd4_df = read_data()
    brd4_df = clean_bb1(brd4_df)
    descriptors_list = get_descriptors()

    descriptor_df = make_descriptor_df(brd4_df, 1)
    for bb in [2, 3]:
        temp_df = make_descriptor_df(brd4_df, bb)
        descriptor_df = pd.concat([descriptor_df, temp_df])
