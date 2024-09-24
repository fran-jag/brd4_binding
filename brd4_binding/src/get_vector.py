



# final_descriptors_list = list(final_df.columns)

# def get_3_tuple(id, column):
#     df_loc_bb1 =  brd4_df.loc[id, 'buildingblock1_smiles']
#     bb1_value = descriptors_bb1_df.loc[df_loc_bb1, column]
    
#     df_loc_bb2 =  brd4_df.loc[id, 'buildingblock2_smiles']
#     bb2_value = descriptors_bb2_df.loc[df_loc_bb2, column]
    
#     df_loc_bb3 =  brd4_df.loc[id, 'buildingblock3_smiles']
#     bb3_value = descriptors_bb3_df.loc[df_loc_bb3, column]

#     return (bb1_value, bb2_value, bb3_value)

# def get_vector_from_id(id):
#     return {column:get_3_tuple(id, column) for column in final_descriptors_list}