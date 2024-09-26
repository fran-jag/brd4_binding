"""
Module to create vector for each id of 3-tuples for each descriptor.

Module that returns a vector of shape (1,210,3) for each id in the
original DataFrame. It creates a 3-tuple for each 2D Descriptor from
RDKit and returns it as a dictionary {descriptor: (tuple)}

"""


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
                       descriptors_list) -> dict:

    return {column: _get_3_tuple(original_df,
                                 id,
                                 descriptor_df,
                                 column)
            for column
            in descriptors_list}


if __name__ == "__main__":
    get_vector_from_id(281065764)
