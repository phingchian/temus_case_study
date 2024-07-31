import os
import argparse
import pandas as pd
from pathlib import Path


def split_ethnicity(pop_table):
    """Splits a population table by ethnicity and saves each group to a separate file.

    Args:
        pop_table (str): Path to the population table file.
    """
    
    df = pd.read_csv(pop_table, sep='\\s+')
    eth_categories = df['ETHNIC_GROUP'].unique()

    for cat in eth_categories:
        cat_df = df[df['ETHNIC_GROUP'] == cat]
        output_path = Path(f'eth_{cat}.txt')
        cat_df.to_csv(output_path, index=False, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split dataframe by ethnicity")
    parser.add_argument('input_file', type=str, help="Path to the population .txt file")
    args = parser.parse_args()

    split_ethnicity(args.input_file)