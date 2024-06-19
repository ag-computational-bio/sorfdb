#!/usr/bin/env python3

import csv
import json
import tarfile
from time import time
from io import BytesIO
from pathlib import Path
from typing import Union, Optional
from argparse import ArgumentParser

import polars as pl
from xopen import xopen

import constants as const


def sequences_of_clusters(cluster_tsv) -> dict[str, list[str]]:
    cluster_sequences: dict[str, list[str]] = dict()
    for protein, cluster in pl.read_csv(cluster_tsv, separator='\t').iter_rows():
        cluster_sequences.setdefault(cluster, [])
        cluster_sequences[cluster].append(protein)
    return cluster_sequences


def extract_description_from_hmms(hmm: str) -> dict[str, str]:
    description_by_id: dict[str, str] = dict()
    hmm_id:  str = ''
    with xopen(hmm, mode='r') as hmm_file:
        for line in hmm_file:
            if line.startswith('//'):
                description_by_id.setdefault(hmm_id, description)
            elif line.startswith('NAME'):
                hmm_id = line.lstrip('NAME').strip()
            elif line.startswith('DESC'):
                description = line.lstrip('DESC').strip()
    return description_by_id


def datapoint_extract_taxonomy(datapoint: dict) -> dict[str, str]:
    required_entries: tuple[str, ...] = ('cluster', 'protein', 'phylum', 'class', 'order', 'family', 'genus', 'species')
    return {k: datapoint[k] for k in required_entries}


def read_cluster_datapoints(fh: Path) -> pl.DataFrame:
    sorfdb_entries: list[dict] = []
    with tarfile.open(fh, "r:gz") as tar_in:
        for entry in tar_in:
            if entry.isfile():
                fh = tar_in.extractfile(entry)
                datapoint: dict = json.load(fh)
                sorfdb_entries.append(datapoint_extract_taxonomy(datapoint))
    return pl.from_dicts(sorfdb_entries).unique()


def read_alignment(cluster: str, alignment_path: Path) -> str:
    file_path: Path = alignment_path.joinpath(f'{cluster}.fasta')
    with open(file_path, mode='r') as fh:
        alignment = fh.read()
    return alignment


def compute_taxonomic_distribution_as_tree(cluster_id: str,
                                           cluster_taxonomy: pl.DataFrame) -> [dict[str, Union[str, int, None]]]:
    """
    Compute the taxonomic distribution of the given small protein family. For each member store only the reported
    taxonomy to store the distribution of the small protein and not the number of annotations.
    :param cluster_id: Cluster ID
    :param cluster_taxonomy: Dataframe with complete cluster taxonomy
    :return: taxonomy as FlatTree
    """
    taxonomy_levels: tuple[str, ...] = ('phylum', 'class', 'order', 'family', 'genus', 'species')
    taxonomy_level_entries: dict[str, dict[str, dict[str, Union[str, int, None]]]] = dict()
    taxonomy_flat_tree: [dict[str, Union[str, int, None]]] = []

    next_id: int = 1
    for sorfdb_datapoint in cluster_taxonomy.filter(pl.col('cluster') == cluster_id).to_dicts():
        parent: Optional[int] = None
        for j, level in enumerate(taxonomy_levels):
            if len(sorfdb_datapoint[level]) > 0:
                taxonomy_level_entries.setdefault(level, dict())
                taxonomy_level_entries[level].setdefault(
                    sorfdb_datapoint[level],
                    {
                        'id': next_id,
                        'rank': level,
                        'label': sorfdb_datapoint[level],
                        'value': 0,
                        'parent': parent
                    }
                )
                taxonomy_level_entries[level][sorfdb_datapoint[level]]['value'] += 1
                parent = taxonomy_level_entries[level][sorfdb_datapoint[level]]['id']
                if taxonomy_level_entries[level][sorfdb_datapoint[level]]['value'] == 1:
                    next_id += 1
            else:
                break

    for level in taxonomy_levels:
        if level in taxonomy_level_entries:
            for tree_node in taxonomy_level_entries[level].values():
                taxonomy_flat_tree.append(tree_node)
    return taxonomy_flat_tree


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Export all sORFs to JSON datapoints for Elasticsearch and to TSV.')
    parser.add_argument('--archive', '-a', type=Path,
                        help='Input sORFdb withClusters as compressed tar archive')
    parser.add_argument('--clusters', '-c', type=Path,
                        help='Input sORFdb small protein clusters as compressed TSV')
    parser.add_argument('--hmm', '-m', type=Path,
                        help='Input sORFdb HMMs as compressed HMM file.')
    parser.add_argument('--alignments', '-l', type=Path,
                        help='Input sORFdb small protein family alignment directory.')
    parser.add_argument('--output', '-o', type=Path, default=Path('./'),
                        help='Output sORFdb tar archive with added cluster information.')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    print('Parse clusters...')
    cluster_df: pl.DataFrame = pl.read_csv(
        args.clusters,
        separator='\t'
    ).groupby('cluster').agg([                                              # 0 cluster id
            pl.col('protein').count().alias('count'),                       # 1 sequenceCount
            pl.col('protein').str.lengths().mean().alias('mean').round(2),  # 2 averageSequenceLength
            pl.col('protein').str.lengths().median().alias('median'),       # 3 medianSequenceLength
        ])

    cluster_sequences: dict[str, list[str]] = sequences_of_clusters(args.clusters)

    print('Parse family descriptions...')
    description_by_id: dict[str, str] = extract_description_from_hmms(args.hmm)

    print('Parse sORFdb entries...')
    cluster_taxonomy: pl.DataFrame = read_cluster_datapoints(args.archive)

    print('Calculate and write cluster statistic...')
    mod_time: float = time()
    with (tarfile.open(args.output.joinpath(f'sorfdb.{const.SORFDB_VER}.clusterStatistic.tar.gz'), 'w:gz') as tar_out,
          xopen(args.output.joinpath(f'sorfdb.{const.SORFDB_VER}.families.tsv.gz'), 'w') as tsv_handle):
        tsv = csv.writer(tsv_handle, delimiter='\t', lineterminator='\n')
        tsv.writerow(('protein', 'family', 'function', 'memberCount', 'averageSequenceLength'))

        for cluster_stats in cluster_df.iter_rows():
            cluster_datapoint: dict = {
                'id': cluster_stats[0],
                'statistics': {
                    'sequenceCount': cluster_stats[1],
                    'averageSequenceLength': cluster_stats[2],
                    'medianSequenceLength': cluster_stats[3],
                    'taxonomy': compute_taxonomic_distribution_as_tree(cluster_stats[0], cluster_taxonomy)
                },
                'function': description_by_id[cluster_stats[0]],
                'alignment': read_alignment(cluster_stats[0], args.alignments)
            }
            datapoint_string = json.dumps(cluster_datapoint, indent=0)
            datapoint_io = BytesIO(bytes(datapoint_string, 'utf-8'))
            info = tarfile.TarInfo(name=f'clusterStats/{cluster_stats[0]}.json')
            info.mtime = mod_time
            info.mode = 0o444
            info.size = len(datapoint_string)
            tar_out.addfile(tarinfo=info, fileobj=datapoint_io)

            tsv.writerows(
                ((protein,
                 cluster_stats[0],
                 description_by_id[cluster_stats[0]],
                 cluster_stats[1],
                 cluster_stats[2]) for protein in cluster_sequences[cluster_stats[0]]))


if __name__ == '__main__':
    main()

# EOF
