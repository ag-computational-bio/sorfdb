#!/usr/bin/env python3

import os
from pathlib import Path
from argparse import ArgumentParser

import polars as pl
# from xopen import xopen


def sequential_melt(df: pl.DataFrame, minsrv: float = 0.0) -> pl.DataFrame:
    """
    Melt the wide Dataframe for sequential for each query to reduce memory usage.
    :param df: Wide Dataframe matrix with subject, query1, query2, ... columns. Values have to SRVs.
    :param minsrv: Report only blast hits with a minimum SRV of minsrv.
    :return: Dataframe in the abc format
    """
    # Create initial dataframe with first query
    abc_df: pl.DataFrame = df.lazy().select(
        pl.col('subject'),
        pl.col(df.columns[1])
    ).melt(
        id_vars='subject',
        variable_name='query',
        value_name='srv'
    ).drop_nulls().filter(
        pl.col('srv') >= minsrv
    ).select(
        [pl.col('query'),
         pl.col('subject'),
         pl.col('srv')]
    ).collect()

    # Create single query melt plans and collect all of them
    lazy_frames = []
    for column in df.columns[2:]:
        lazy_frames.append(
            df.lazy().select(
                pl.col('subject'),
                pl.col(column)
            ).melt(
                id_vars='subject',
                variable_name='query',
                value_name='srv'
            ).drop_nulls().filter(
                pl.col('srv') >= minsrv
            ).select(
                [pl.col('query'),
                 pl.col('subject'),
                 pl.col('srv')]
            )
        )

    for tmp_df in pl.collect_all(lazy_frames):
        abc_df.vstack(tmp_df, in_place=True)
    return abc_df.rechunk()


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Import the diamond  the ')
    parser.add_argument('--input', '-i', type=Path, help='Input a BLAST TSV file with outfmt 6')
    parser.add_argument('--output', '-o', type=Path, default=Path('./'), help='Output folder path')
    parser.add_argument('--mode', '-m', default='blast', type=str, help='Mode for "diamond" or "blast"')
    parser.add_argument('--minsrv', '-n', type=int, default=0,
                        help='Report only blast hits with a minimum SRV of value 30 (= 3.0 * 10).')
    parser.add_argument('--id', '-d', type=int, default=0,
                        help='Report only alignments above the given percentage of sequence identity.')
    parser.add_argument('--query_cover', '-q', type=int, default=0,
                        help='Report only alignments above the given percentage of query cover.')
    parser.add_argument('--subject_cover', '-s', type=int, default=0,
                        help='Report only alignments above the given percentage of subject cover.')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    os.makedirs(args.output, exist_ok=True)

    if args.mode not in {'diamond', 'blast'}:
        raise Exception('Supported modes: "diamond" or "blast"')

    print(f'Settings:\n'
          f'\tmode: {args.mode}\n'
          f'\tid: {args.id}\n'
          f'\tquery_cover: {args.query_cover}\n'
          f'\tsubject_cover: {args.subject_cover}\n'
          f'\tminsrv: {args.minsrv}\n')

    print(f'Load {args.mode} results...')
    if args.mode == 'diamond':
        df = pl.read_csv(
            args.input,
            has_header=False,
            new_columns=['query', 'subject', 'bitscore', 'evalue'],
            separator='\t'
        ).drop('evalue').groupby(
            'query',
            'subject'
        ).agg(pl.col('bitscore').max())
    else:  # blast
        df = pl.scan_csv(
            args.input,
            has_header=False,
            new_columns=['query', 'subject', 'identity', 'alignment_length', 'bitscore', 'evalue'],
            dtypes={'bitscore': pl.Float64},
            separator='\t'
        ).filter(
            (pl.col('identity') >= args.id) &
            (pl.col('alignment_length') / pl.col('query').str.lengths() * 100 >= args.query_cover) &
            (pl.col('alignment_length') / pl.col('subject').str.lengths() * 100 >= args.subject_cover)
        ).select(
            pl.col('query'),
            pl.col('subject'),
            pl.col('bitscore')
        ).groupby(
            'query',
            'subject'
        ).agg(pl.col('bitscore').max()).collect()

    print('Pruning...')  # exclude self-hit only
    df = df.filter(
        ~pl.col('query').is_in(
            df.groupby('query').count().filter(pl.col('count') == 1).filter(
                # unique query is in unique subject == self-hit
                pl.col('query').is_in(
                    df.groupby('subject').count().filter(pl.col('count') == 1).select(pl.col('subject')).to_series())
            ).select('query').to_series()
        )
    ).rechunk()

    print('Transform to matrix...')
    df = df.pivot(index='subject',
                  columns='query',
                  values='bitscore')  # .fill_null(0) do not fill nulls to create a sparse matrix

    # print('Write bitscore matrix...')
    # with xopen(args.output.joinpath('bitscore.matrix.tsv.gz'), 'w', compresslevel=9, threads=args.threads) as of:
    #     of.write(df.write_csv(separator='\t'))

    print('Calculate and rescale SRVs...')
    df = df.lazy().select(
        pl.col('subject'),
        *(pl.col(column) / pl.col(column).max() * 100 for column in df.columns[1:])  # max per column == query self-hit
    ).collect()

    # print('Write SRV matrix...')
    # with xopen(args.output.joinpath('srv.matrix.tsv.gz'), 'w', compresslevel=9, threads=args.threads) as of:
    #     of.write(df.write_csv(separator='\t'))

    print('Reformat to abc format...')
    df = sequential_melt(df, minsrv=float(args.minsrv))

    # Uses too much RAM
    # df = df.lazy().melt(
    #     id_vars='subject',
    #     variable_name='query',
    #     value_name='srv'
    # ).drop_nulls().filter(
    #     pl.col('srv') >= args.minsrv
    # ).select(
    #     [pl.col('query'),
    #      pl.col('subject'),
    #      pl.col('srv')]
    # ).collect()  # .sort(['query', 'subject'])

    print('Pruning & rounding...')  # exclude all self hits
    df = df.filter(
        pl.col('query') != pl.col('subject')
    ).select(
        pl.col('query'),
        pl.col('subject'),
        pl.col('srv').round(4)
    )

    print('Write SRV.abc...')
    # with xopen(args.output.joinpath('srv.tsv.gz'), 'w', compresslevel=9, threads=args.threads) as of:
    #     of.write(df.write_csv(separator='\t'))

    df.write_csv(args.output.joinpath('srv.tsv'), separator='\t', has_header=False)


if __name__ == '__main__':
    main()

# EOF
