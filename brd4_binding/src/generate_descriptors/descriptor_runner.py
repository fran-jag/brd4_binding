import pandas as pd

from rdkit.rdBase import BlockLogs

from pipeline.collection import read_data
from pipeline.create_descriptor_df import make_descriptor_df, get_descriptors
from pipeline.create_3tuple_vector import get_vector_from_id
from pipeline.remove_protecting_group import clean_bb1


def main():
    # Silence RDKit deprecation warning
    _block = BlockLogs()

    brd4_df = read_data("/home/papafrita/Projects/brd4_binding/data/brd4_5p.parquet")
    brd4_df = clean_bb1(brd4_df)
    descriptors_list = get_descriptors()

    descriptor_df = make_descriptor_df(brd4_df, 1)
    for bb in [2, 3]:
        temp_df = make_descriptor_df(brd4_df, bb)
        descriptor_df = pd.concat([descriptor_df, temp_df],
                                  axis=0)
    descriptor_df.drop_duplicates(inplace=True)
    ids_test = brd4_df.index
    output_df = pd.DataFrame(columns=descriptors_list)
    id_counter = 0
    temp_dict = {}

    for id in ids_test:
        temp_row = get_vector_from_id(original_df=brd4_df,
                                      id=id,
                                      descriptor_df=descriptor_df,
                                      descriptors_list=descriptors_list,
                                      )
        temp_dict[id] = temp_row

        if (id_counter % 2500) == 0:
            print("{} of {} processed.".format(id_counter, brd4_df.shape[0]))
        id_counter += 1

    output_df = pd.DataFrame.from_dict(temp_dict, orient="index")

    return output_df.join(brd4_df['binds'])


if __name__ == "__main__":
    output_df = main()
    output_df.to_parquet("/home/papafrita/Projects/brd4_binding/data/output_5p.parquet")
