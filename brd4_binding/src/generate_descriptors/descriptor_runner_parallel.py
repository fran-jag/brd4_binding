import pandas as pd

from rdkit.rdBase import BlockLogs

from pipeline.collection import read_data
from pipeline.create_descriptor_df import make_descriptor_df, get_descriptors
from pipeline.create_3tuple_vector import get_vector_from_id
from pipeline.remove_protecting_group import clean_bb1

from multiprocessing import Pool, cpu_count
import numpy as np


def process_chunk(ids,
                  original_df,
                  descriptor_df,
                  descriptor_list,
                  ):
    chunk_result = {}
    for id in ids:
        chunk_result[id] = get_vector_from_id(original_df,
                                              id,
                                              descriptor_df,
                                              descriptor_list,
                                              )
    return chunk_result


def parallelize_processing(ids, func, n_cores=None):
    if n_cores is None:
        n_cores = cpu_count()

    ids_split = np.array_split(ids, n_cores)
    pool = Pool(n_cores)
    results = pool.map(func, ids_split)
    pool.close()
    pool.join()

    combined_result = {}
    for result in results:
        combined_result.update(result)
    return combined_result


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

    # Parallelize the processing
    all_ids = brd4_df.index[:50000]
    parallel_result = parallelize_processing(all_ids,
                                             process_chunk(original_df=brd4_df,
                                                           descriptor_df=descriptor_df,
                                                           descriptor_list=descriptors_list,
                                                           )
                                             )

    output_df = pd.DataFrame.from_dict(parallel_result, orient='index')
    output_df.to_pickle("/home/papafrita/Projects/brd4_binding/data/output_parallel_50k")


if __name__ == "__main__":
    output_df = main()
