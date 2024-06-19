#!/usr/bin/env python3

from pathlib import Path
from argparse import ArgumentParser
from itertools import product

import polars as pl

import constants as const


def read_mcl_cls(file: Path) -> dict:
    cluster_map: dict[str, int] = dict()
    with open(file, 'r') as infile:
        for i, line in enumerate(infile):
            proteins: list[str] = line.strip().split()
            for protein, cls in product(proteins, [str(i)]):
                cluster_map[protein] = cls
    return cluster_map


def load_ribosomal_proteins_from_db(db_path: Path) -> pl.DataFrame:
    """
    Extract protein sequences and products for ribosomal proteins.
    :param db_path: path to sorfdb tsv
    :return:
    """
    if str(db_path).endswith('.gz'):
        opener = pl.read_csv
    else:
        opener = pl.scan_csv

    db_df = opener(
        db_path,
        separator='\t'
    ).select(
        pl.col('protein'),
        pl.col('product')
    ).filter(
        (pl.col('product').str.contains('ibosomal')) &
        (pl.col('protein').str.lengths() <= const.HMM_MAX_SPROT_LENGTH) &
        (~pl.col('product').str.to_lowercase().is_in(pl.Series('hypothetical_product_names',
                                                     [term for term in const.HYPOTHETICAL_PRODUCT_NAMES])))
    ).with_columns(
        pl.col('product').str.to_lowercase()
    ).with_columns(
        pl.col('product').str.extract(
            r"(\d{1,2}\w{1}|\wsu|small subunit).{1}ribosomal.{1}(subunit.{1})?protein.{1}(subunit.{1})?(\w{1}\d{1,2}|\w{3})",
            0
        ).str.replace(r"\+|%",
                      ' ').str.replace(r"\+|%",
                                       ' ').str.replace(r"\+|%",
                                                        ' ').str.replace(r"ssu",
                                                                         '30s').str.replace(r"lsu",
                                                                                            '50s').str.replace('subunit ', '')
    ).with_columns(
        pl.when(
            pl.col('product').str.lengths() > 0
        ).then(
            pl.lit(1)
        ).otherwise(
            pl.lit(0)
        ).alias('ribosomal')
    )

    if not str(db_path).endswith('.gz'):
        db_df = db_df.collect()
    return db_df.rename({'protein': 'query'})


def export_ribosomal_clusters(graph_df: pl.DataFrame, inflations: list[str], outpath: Path) -> pl.DataFrame:
    """
    Calculate the ratio of clusters which cover a family of ribosomal proteins.
    :param graph_df: graph Dataframe in abc format
    :param inflations: clusters for different inflations
    :param outpath: path to the output file
    :return: Dataframe
    """
    coverage_df: pl.DataFrame = graph_df.lazy().filter(
            pl.col('ribosomal') == 1
        ).groupby(
            inflations[0], 'product'
        ).count().groupby('product').agg(
            (pl.col('count') / pl.col('count').sum() * 100).round(2).sort().alias(f'coverage-{inflations[0]}')
        ).with_columns(
            pl.col(f'coverage-{inflations[0]}').cast(pl.List(pl.Utf8)).list.join(", ")
        ).sort('product').collect()

    lazy_frames: list[pl.LazyFrame] = []
    for inflation in inflations[1:]:
        lazy_frames.append(
            graph_df.lazy().filter(
                pl.col('ribosomal') == 1
            ).groupby(
                inflation, 'product'
            ).count().groupby('product').agg(
                (pl.col('count') / pl.col('count').sum() * 100).round(2).sort().alias(f'coverage-{inflation}')
            ).with_columns(
                pl.col(f'coverage-{inflation}').cast(pl.List(pl.Utf8)).list.join(", ")
            ).sort('product')
        )

    for tmp_df in pl.collect_all(lazy_frames):
        coverage_df = coverage_df.join(tmp_df, on='product', how='left')
    coverage_df.write_csv(outpath.joinpath('ribosomal_cluster_coverage.tsv'), separator='\t')
    return coverage_df.rechunk()


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Import a mcl abc file and tab files with clusters and export them for use in cytoscape.')
    parser.add_argument('--abc', '-a', type=Path, help='Input mcl tab file with clusters')
    parser.add_argument('--sorfdb', '-s', type=Path, help='sORFDB TSV file')
    parser.add_argument('--autoclust', '-u', type=Path, help='MCL clusters for inflation values chosen for binned subgraphs.')
    parser.add_argument('--cluster', '-c', nargs='+', type=list[list[str]], dest='cluster',
                        help='MCL clusters for different inflation values.')
    parser.add_argument('--output', '-o', type=Path, default=Path('./'), help='Output folder path')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    graph_df = pl.read_csv(
        args.abc,
        has_header=False,
        new_columns=['query', 'subject', 'srv'],
        separator='\t'
    )

    inflations: list[str] = ['q-autoclust']
    clusters = read_mcl_cls(args.autoclust)
    graph_df = graph_df.with_columns(
        pl.col('query').map_dict(clusters).alias('q-autoclust')
    )

    for clustering in args.cluster:
        clustering: Path = Path(''.join(clustering))
        inflation: str = '.'.join(str(clustering).split('/')[-1].split('.')[-3:-1])
        inflations.append(f'q-{inflation}')
        clusters = read_mcl_cls(clustering)

        graph_df = graph_df.with_columns(
            pl.col('query').map_dict(clusters).alias(f'q-{inflation}')
        )

    # export_ribosomal_clusters(
    #     graph_df.join(sorfdb_df, on='query', how='left'),
    #     inflations,
    #     args.output
    # )

    graph_df.write_csv(args.output.joinpath('srv.clustered.tsv'),
                       separator='\t')


if __name__ == '__main__':
    main()

# EOF
