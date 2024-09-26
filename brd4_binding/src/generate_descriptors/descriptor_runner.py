import pandas as pd
import pickle
import time

from rdkit.rdBase import BlockLogs

from pipeline.collection import read_data
from pipeline.create_descriptor_df import make_descriptor_df, get_descriptors
from pipeline.create_3tuple_vector import get_vector_from_id
from pipeline.remove_protecting_group import clean_bb1


def main():
    # Silence RDKit deprecation warning
    _block = BlockLogs()

    brd4_df = read_data("/home/papafrita/Projects/" +
                        "brd4_binding/data/brd4_5p.parquet")
    brd4_df = clean_bb1(brd4_df)
    descriptors_list = get_descriptors()

    descriptor_df = make_descriptor_df(brd4_df)

    ids_test = brd4_df.index[:100000]

    start_time = time.time()
    print("Starting data processing...")
    print(time.strftime("Start time: %H:%M:%S",
                        time.localtime(start_time)))

    id_counter = 0
    temp_dict = {}
    failed_ids = []
    for id in ids_test:
        try:
            temp_row = get_vector_from_id(original_df=brd4_df,
                                          id=id,
                                          descriptor_df=descriptor_df,
                                          descriptors_list=descriptors_list,
                                          )
            temp_dict[id] = temp_row
        except KeyError:
            failed_ids.append(id)

        if (id_counter % 5000) == 0:

            current_time = time.time()
            delta_time = current_time - start_time

            print("{} of {} processed.".format(id_counter, len(ids_test)))
            print("{:.1f}%".format(id_counter/len(ids_test)*100))

            if id_counter > 0:
                factor = len(ids_test)/id_counter
                print(time.strftime("ETA: %M minutes",
                                    time.localtime(delta_time*factor)))
            print("-"*10)

        id_counter += 1

    output_df = pd.DataFrame.from_dict(temp_dict, orient="index")

    end_time = time.time()
    print("End of data processing. Saving...")
    print(time.strftime("End time: %H:%M:%S",
                        time.localtime(end_time)))

    if len(failed_ids) > 0:
        print("{} IDs failed".format(len(failed_ids)))
        print("\nSaving list in /data/failed_ids.pickle")

        with open("/home/papafrita/Projects/" +
                  "brd4_binding/data/failed_ids.pickle",
                  "wb") as file:
            pickle.dump(failed_ids, file)
        print("DataFrame saved.")
    else:
        print("No IDs failed. Saving DataFrame.")

    return output_df.join(brd4_df['binds'])


if __name__ == "__main__":
    output_df = main()
    output_df.to_parquet("/home/papafrita/Projects/" +
                         "brd4_binding/data/output_100k.parquet")
