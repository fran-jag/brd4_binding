"""
Module to remove any protecting group form a buildingblock.

Module that removes a specificied protecting group from a
buildingblock in the original data. Returns a pandas.DataFrame
"""

import pandas as pd


def clean_bb1(original_data: pd.DataFrame,
              buildingblock=1,
              prot_group='CC1c2ccccc2-c2ccccc21',
              ) -> pd.DataFrame:

    buildingblock_col = "buildingblock{}_smiles".format(buildingblock)

    for index in original_data.index:
        old_str = original_data.loc[index, buildingblock_col]
        new_val = old_str.replace(prot_group, '')
        original_data.at[index, buildingblock_col] = new_val

    return original_data
