#!/usr/bin/env python3

import sys
import click
import numpy as np
import pandas as pd
from pathlib import Path


def nan_divide(x, y):
    try:
        return x / y
    except:
        return np.nan


@click.command()
@click.argument('grouped-samples', type=click.Path(exists=True))
def main(grouped_samples):
    """ Calculate the sharing percentages between tumor and spleen samples
    of the same mouse

    Input: 5 column csv: mouse_id,spleen_sample,spleen_path,tumor_sample,tumor_path,
    where spleen_path and tumor point to the CDR3 peptide files.

    See assets/grouped_samples.csv
    """
    
    df = pd.read_csv(grouped_samples)

    new_rows = []
    for _, row in df.iterrows():

        df_tumor = pd.read_csv(row['tumor_path'], header=None, names=['cdr3', 'Tumor'], index_col='cdr3')
        df_spleen = pd.read_csv(row['spleen_path'], header=None, names=['cdr3', 'Spleen'], index_col='cdr3')
        df_merged = pd.merge(df_tumor, df_spleen, left_index=True, right_index=True, how='inner')

        new_row = {}
        new_row['chain'] = row['chain']
        new_row['mouse_id'] = row['mouse_id']
        new_row['Tumor Sample'] = row['tumor_sample']
        new_row['Tumor uCDR3'] = len(df_tumor)
        new_row['Spleen Sample'] = row['spleen_sample']
        new_row['Spleen uCDR3'] = len(df_spleen)
        new_row['Shared uCDR3'] = len(df_merged)
        new_row['Fraction of CDR3s from Spleen in Tumor'] = nan_divide(len(df_merged), len(df_spleen))
        new_row['Fraction of CDR3s from Tumor in Spleen'] = nan_divide(len(df_merged), len(df_tumor))
        new_row['Fraction of Expression of Spleen in Tumor'] = nan_divide(df_merged['Spleen'].sum(), df_spleen['Spleen'].sum())
        new_row['Fraction of Expression of Tumor in Spleen'] = nan_divide(df_merged['Tumor'].sum(), df_tumor['Tumor'].sum())
        new_rows.append(new_row)

    df_final = pd.DataFrame(new_rows)
    df_final.to_csv(sys.stdout, index=False)


if __name__ == "__main__":
    main()
