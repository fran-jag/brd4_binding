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


def read_data(path="../../data/brd4_5p.parquet",
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
    print("Reading parquet file...")
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
def get_descriptors() -> list:
    descriptors_list = [x[0] for x
                        in Descriptors._descList
                        ]
    return descriptors_list


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


def _get_3_tuple(original_df,
                 id,
                 descriptor_df,
                 column,
                 ) -> tuple:
    df_loc_bb1 = original_df.loc[id, 'buildingblock1_smiles']
    bb1_value = descriptor_df.loc[df_loc_bb1, column]

    df_loc_bb2 = original_df.loc[id, 'buildingblock2_smiles']
    bb2_value = descriptor_df.loc[df_loc_bb2, column]

    df_loc_bb3 = original_df.loc[id, 'buildingblock3_smiles']
    bb3_value = descriptor_df.loc[df_loc_bb3, column]

    return (bb1_value, bb2_value, bb3_value)


def get_vector_from_id(original_df,
                       id,
                       descriptor_df,
                       descriptors_list):
    return {column: _get_3_tuple(original_df,
                                 id,
                                 descriptor_df,
                                 column)
            for column
            in descriptors_list}


if __name__ == "__main__":
    # Silence RDKit deprecation warning
    block = BlockLogs()

    brd4_df = read_data()
    brd4_df = clean_bb1(brd4_df)
    descriptors_list = get_descriptors()

    descriptor_df = make_descriptor_df(brd4_df, 1)
    for bb in [2, 3]:
        temp_df = make_descriptor_df(brd4_df, bb)
        descriptor_df = pd.concat([descriptor_df, temp_df],
                                  axis=0)
    descriptor_df.drop_duplicates(inplace=True)

    ids_test = brd4_df.index[:100]
    output_df = pd.DataFrame(columns=descriptors_list,
                             index='id')
    for id in ids_test:
        temp_row = get_vector_from_id(original_df=brd4_df,
                                      id=id,
                                      descriptor_df=descriptor_df,
                                      descriptors_list=descriptors_list,
                                      )
        output_df.loc[id] = temp_row
    print(output_df.head())
