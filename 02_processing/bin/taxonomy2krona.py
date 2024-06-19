#!/usr/bin/env python3

from pathlib import Path
from argparse import ArgumentParser

import polars as pl


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Import the taxonomy of sORFdb as TSV with .')
    parser.add_argument('--input', '-i', type=Path, help='Input path to the sORFdb file.')
    parser.add_argument('--output', '-o', type=Path, default=Path('./'), help='Output folder path')
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()
    column_names: list[str] = ['phylum', 'class', 'order', 'family', 'genus', 'protein']

    # total
    pl.scan_csv(
        args.input,
        has_header=False,
        new_columns=column_names,
        separator='\t'
    ).select(
        *[pl.col(level) for level in column_names[:-1]]
    ).filter(
        ~pl.col('phylum').str.starts_with('Candidatus') &
        ~pl.col('phylum').str.starts_with('Candidate')
    ).groupby(
        'phylum', 'class', 'order', 'family', 'genus'
    ).count().with_columns(
        pl.lit('Bacteria').alias('kingdom')
    ).select(
        pl.col('count'),
        pl.col('kingdom'),
        *[pl.col(colname) for colname in column_names[:-1]]
    ).sink_csv(
        args.output.joinpath('taxonomy.total.tsv'),
        has_header=False,
        separator='\t',
    )

    # non-redundant
    pl.scan_csv(
        args.input,
        has_header=False,
        new_columns=column_names,
        separator='\t'
    ).filter(
        ~pl.col('phylum').str.starts_with('Candidatus') &
        ~pl.col('phylum').str.starts_with('Candidate')
    ).unique().groupby(
        'phylum', 'class', 'order', 'family', 'genus'
    ).count().with_columns(
        pl.lit('Bacteria').alias('kingdom')
    ).select(
        pl.col('count'),
        pl.col('kingdom'),
        *[pl.col(colname) for colname in column_names[:-1]]
    ).sink_csv(
        args.output.joinpath('taxonomy.non-redundant.tsv'),
        has_header=False,
        separator='\t',
    )


if __name__ == '__main__':
    main()

# EOF
