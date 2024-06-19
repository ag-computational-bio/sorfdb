#!/usr/bin/env python3

import os
import csv
from pathlib import Path
from argparse import ArgumentParser

import polars as pl
from xopen import xopen


def load_query_and_subject_seqs(abc_file: Path) -> tuple[set[str], set[str], set[str], set[str]]:
    """
    Extract all as well as unique query and subject sequences from a abc TSV file.
    a = query, b = subject, c = value
    :param abc_file: path to abc TSV file
    :return: unique queries, unique subjects
    """
    query_counter: dict[str, int] = dict()
    subject_counter: dict[str, int] = dict()

    with xopen(abc_file, 'r') as file_in:
        tsvreader = csv.reader(file_in, delimiter='\t')
        for hit in tsvreader:
            query, subject, _ = hit
            query_counter.setdefault(query, 0)
            query_counter[query] += 1
            subject_counter.setdefault(subject, 0)
            subject_counter[subject] += 1

    unique_queries: set[str] = set(seq for seq, count in query_counter.items() if count == 1)
    unique_subjects: set[str] = set(seq for seq, count in subject_counter.items() if count == 1)
    return unique_queries, unique_subjects, set(query_counter.keys()), set(subject_counter.keys())


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Import a SRV abc file and exclude single hits.')
    parser.add_argument('--input', '-i', type=Path, help='Input SRV values in abc format')
    parser.add_argument('--output', '-o', type=Path, default=Path('./'), help='Output folder path')
    parser.add_argument('--low-memory', '-l', action='store_true', dest='low_memory',
                        help='Low memory option for single hit pruning.')
    parser.add_argument('--compress', '-c', action='store_true',
                        help='Compress the output.')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    os.makedirs(args.output, exist_ok=True)

    if args.low_memory:
        unique_queries, unique_subjects, queries, subjects = load_query_and_subject_seqs(args.input)
        queries = queries - unique_queries
        subjects = subjects - unique_subjects

        if args.compress:
            outpath: Path = args.output.joinpath('pruned.srv.tsv.gz')
        else:
            outpath: Path = args.output.joinpath('pruned.srv.tsv')

        counter: int = 0
        with xopen(args.input, 'r') as file_in, xopen(outpath, 'w', compresslevel=9, threads=args.threads) as file_out:
            tsvreader = csv.reader(file_in, delimiter='\t')
            tsvwriter = csv.writer(file_out, delimiter='\t', lineterminator='\n')
            for hit in tsvreader:
                # exclude single hits and exclude non-bidirectional hits (query has to be in subjects and vice versa)
                if hit[0] in unique_queries or hit[1] in unique_subjects:
                    counter += 1
                elif hit[0] in subjects and hit[1] in queries:
                    tsvwriter.writerow(hit)
                else:
                    counter += 1

                # if (hit[0] not in unique_queries and hit[1] not in unique_subjects) and \
                #         (hit[0] in subjects and hit[1] in queries):

        print('Pruned:', counter)
        del unique_queries, unique_subjects, counter
    else:
        # unique_queries: pl.Series = pl.scan_csv(
        #     args.input,
        #     has_header=False,
        #     new_columns=['query', 'subject', 'srv'],
        #     dtypes={'srv': pl.Float64},
        #     low_memory=True,
        #     separator='\t'
        # ).select(pl.col('query')).groupby('query').count().filter(pl.col('count') == 1).select(
        #     pl.col('query')
        # ).collect().to_series()
        #
        # unique_subjects: pl.Series = pl.scan_csv(
        #     args.input,
        #     has_header=False,
        #     new_columns=['query', 'subject', 'srv'],
        #     dtypes={'srv': pl.Float64},
        #     low_memory=True,
        #     separator='\t'
        # ).select(
        #     pl.col('subject')
        # ).groupby('subject').count().filter(pl.col('count') == 1).select(
        #     pl.col('subject')
        # ).collect().to_series()
        #
        # # drop edges with only one connection
        # df = pl.scan_csv(
        #     args.input,
        #     has_header=False,
        #     new_columns=['query', 'subject', 'srv'],
        #     dtypes={'srv': pl.Float64},
        #     low_memory=True,
        #     separator='\t'
        # ).filter(
        #     ~(pl.col('query').is_in(unique_queries)) &
        #     ~(pl.col('subject').is_in(unique_subjects))
        # ).collect()
        # del unique_queries, unique_subjects
        #
        # unique_queries: set[str] = set(pl.scan_csv(
        #     args.input,
        #     has_header=False,
        #     new_columns=['query', 'subject', 'srv'],
        #     dtypes={'srv': pl.Float64},
        #     low_memory=True,
        #     separator='\t'
        # ).select(pl.col('query')).unique(
        #     keep='none'
        # ).collect().to_series().to_list())
        #
        # unique_subjects: set[str] = set(pl.scan_csv(
        #     args.input,
        #     has_header=False,
        #     new_columns=['query', 'subject', 'srv'],
        #     dtypes={'srv': pl.Float64},
        #     low_memory=True,
        #     separator='\t'
        # ).select(
        #     pl.col('subject')
        # ).unique(
        #     keep='none'
        # ).collect().to_series().to_list())

        df = pl.scan_csv(
            args.input,
            has_header=False,
            new_columns=['query', 'subject', 'srv'],
            dtypes={'srv': pl.Float64},
            low_memory=True,
            separator='\t'
        ).filter(
            ~(pl.col('query').is_unique()) &
            ~(pl.col('subject').is_unique())
        ).collect()

        if args.compress:
            with xopen(args.output.joinpath('pruned.srv.tsv.gz'), 'w', compresslevel=9,
                       threads=args.threads) as file_out:
                file_out.write(df.write_csv(separator='\t', has_header=False))
        else:
            df.write_csv(args.output.joinpath('pruned.srv.tsv'), separator='\t', has_header=False)


if __name__ == '__main__':
    main()

# EOF
