"""
Module for reading data.

Read data from parquet file and returns pandas.DataFrame
with a specified column as index.

Dependencies:
    - pandas

"""
import pandas as pd


def read_data(path="../data/brd4_5p.parquet", index='id') -> pd.DataFrame:
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