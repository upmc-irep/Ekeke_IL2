#!/usr/bin/env python3

import re
import sys
import click
import pandas as pd
import numpy as np
from pathlib import Path
from skbio.diversity.alpha import pielou_e


def d50_head(df: pd.DataFrame) -> float:
    """ What percentage of unique cdr3s of the head make up
    50% of the total reads in the head?

    d50_head = |F| / |H|,
        where H is set of UCDR3s with frequency >= the mean of the frequency distribution
        where F is the set of uCDR3s that make up the top 50% of reads in H
    """

    head = df[df['counts'] >= df['counts'].mean()]
    half_of_head = head[head['counts'].cumsum() <= head['counts'].sum() * 0.50]
    d50 = ((len(half_of_head) + 1) / len(head)) * 100.0

    return d50


@click.command()
@click.argument('files', type=click.Path(exists=True), nargs=-1, required=True)
def main(files):
    """ Calculate diversity metrics from a list of CDR3 peptide files """ 

    rows = []
    for file in files:
        name = re.sub('_*(TRA|TRB|TRD|TRG|IGH|IGK|IGL)', '', Path(file).stem)

        df = pd.read_csv(file, header=None, names=['cdr3', 'counts'])
        df.sort_values('counts', inplace=True, ascending=False)
        df['freq'] = df['counts'] / df['counts'].sum()

        row = {}
        row['sample'] = name
        row['uCDR3'] = len(df)
        row['d50'] = d50_head(df)
        row['1-Pielou Clonality'] = 1 - pielou_e(df['counts'])
        row['Shannon'] = (-(df['freq'] * np.log2(df['freq'])).sum())
        row['Normalized Shannon'] = row['Shannon'] / np.log2(len(df))
        row['Alpha Index'] = (df.loc[0]['counts'] / df['counts'].sum()) * 100.0
        rows.append(row)

    df_final = pd.DataFrame(rows)
    df_final.to_csv(sys.stdout, index=False)


if __name__ == "__main__":
    main()