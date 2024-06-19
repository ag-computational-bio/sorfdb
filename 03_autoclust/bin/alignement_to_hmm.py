#!/usr/bin/env python3

import csv
from math import floor
from pathlib import Path
from typing import Iterator, Optional, Union
from collections import namedtuple
from argparse import ArgumentParser

import pyhmmer
import polars as pl
from xopen import xopen

import extended_io as eio
from constants import HMM_MAX_SPROT_LENGTH, HYPOTHETICAL_PRODUCT_NAMES, SORFDB_VER


def load_db(db_path: Path) -> pl.DataFrame:
    """
    Import the protein sequences and products from the (compressed) sorfdb.
    :param db_path: path to sorfdb file
    :return: Dataframe with proteins and products from sorfdb
    """
    compressed: bool = str(db_path).endswith('.gz')
    opener: Union[pl.read_csv, pl.scan_csv]
    if compressed:
        opener = pl.read_csv
    else:
        opener = pl.scan_csv

    db_df: Union[pl.DataFrame, pl.LazyFrame]
    db_df = opener(
        db_path,
        separator='\t'
    ).select(
        pl.col('protein'),
        pl.col('product')
    ).filter(
        (pl.col('protein').str.lengths() <= HMM_MAX_SPROT_LENGTH) &
        (~pl.col('product').str.to_lowercase().is_in(pl.Series('hypothetical_product_names',
                                                               [term for term in HYPOTHETICAL_PRODUCT_NAMES])))
    ).with_columns(
        pl.col('product').apply(lambda x: x.split(' (')[0])  # Cut off reference ids of revised hypothetical products
    )

    if compressed:
        return db_df
    return db_df.collect()


def create_cluster_name(proteins: pl.Series, db_frame: pl.DataFrame) -> bytes:
    """
    Create a name of the given proteins based on a major vote on their protein product.
    :param proteins: proteins of a cluster
    :param db_frame: database with protein descriptions
    :return:
    """
    # one most common product
    name_frame = db_frame.lazy().filter(
        pl.col('protein').is_in(proteins)
    ).select(
        pl.col('product')
    ).groupby('product').count().filter(
        pl.col('count') == pl.col('count').max()
    ).collect()

    name: bytes
    if name_frame.shape[0] == 1:
        name = bytes(name_frame.to_series().to_list()[0], 'UTF-8')
    else:
        # one most common product in lower case
        name_frame = db_frame.lazy().filter(
            pl.col('protein').is_in(proteins)
        ).select(
            pl.col('product').str.to_lowercase()
        ).groupby('product').count().filter(
            pl.col('count') == pl.col('count').max()
        ).collect()
        if name_frame.shape[0] == 1:
            name = bytes(name_frame.to_series().to_list()[0], 'UTF-8')
        else:
            name = bytes('', 'UTF-8')
    return name.strip()


def msa_proteins_to_digital_seqs(seqs: tuple[bytes, ...],
                                 alphabet: pyhmmer.easel.Alphabet = pyhmmer.easel.Alphabet.amino()) -> list[pyhmmer.easel.DigitalSequence]:
    """
    Convert proteins to a list of pyhmmer DigitalSequences.
    :param seqs: byte strings of protein sequences
    :param alphabet: alphabet of the sequences
    :return: list of pyhmmer DigitalSequences
    """
    sequences: list[pyhmmer.easel.DigitalSequence] = []
    for aa in seqs:
        sequences.append(
            pyhmmer.easel.TextSequence(sequence=aa.decode('UTF-8'),
                                       name=aa).digitize(alphabet)
        )
    return sequences


def fasta_proteins_to_digital_seqs(fasta_file: Path) -> list[pyhmmer.easel.DigitalSequence]:
    """
    Convert protein sequences of a FASTA file to a list of digital pyhmmer Sequences.
    :param fasta_file: path to FASTA file
    :return: list of pyhmmer DigitalSequences from the sORF proteins
    """
    sequences: list[pyhmmer.easel.DigitalSequence] = []
    for _, aa in eio.parse_fasta(fasta_file):
        sequences.append(
            pyhmmer.easel.TextSequence(sequence=aa,
                                       name=aa.encode('UTF-8')).digitize(pyhmmer.easel.Alphabet.amino())
        )
    return sequences


def perform_hmmsearch(proteins: list[pyhmmer.easel.DigitalSequence],
                      hmm: pyhmmer.plan7.HMM,
                      cutoff: Optional[str] = None,
                      threads: int = 1) -> list[namedtuple]:
    """
    Perform a hmmsearch on the given proteins with the HMM profile(s).
    :param proteins: list of proteins in pyhmmer DigitalSequence format
    :param hmm: HMM profile(s)
    :param cutoff: cutoff threshold
    :param threads: number of threads
    :return: list of best hits (namedtuple)
    """
    Result: namedtuple = namedtuple('Result', ['score', 'best_domain_score', 'included', 'sequence'])
    results: list[namedtuple] = []

    search: Iterator[pyhmmer.plan7.TopHits | pyhmmer.plan7.TopHits]
    search = pyhmmer.hmmsearch(hmm, proteins, bit_cutoffs=cutoff, cpus=threads) if cutoff else pyhmmer.hmmsearch(hmm, proteins, cpus=threads)

    for top_hits in search:
        for hit in top_hits:
            if hit.included:
                results.append(Result(score=hit.score,
                                      best_domain_score=hit.best_domain.score,
                                      included=hit.included,
                                      sequence=hit.name))
    return results


def calculate_gathering_cutoff(hits) -> tuple[float, float]:
    """
    Calculate the gathering cutoffs from the lowest scoring included sequence.
    :param hits: hits of the proteins against their HMM
    :return: gathering cutoffs (query, domain)
    """
    try:
        gathering_score = min(hit.score for hit in hits if hit.included)
        gathering_domscore = min(hit.best_domain_score for hit in hits if hit.included)
    except ValueError:
        gathering_score, gathering_domscore = 0.0, 0.0
    return floor(gathering_score), floor(gathering_domscore)


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Import a SRV abc file and exclude single hits.')
    parser.add_argument('--alignments', '-a', type=Path, help='Input path to directory with alignment files *.aln')
    parser.add_argument('--proteins', '-p', type=Path, help='Input complete sORF FASTA file.')
    parser.add_argument('--sorfdb', '-s', type=Path, help='Input sorfdb.tsv(.gz)')
    parser.add_argument('--output', '-o', type=Path, default=Path('./'), help='Output folder path')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    print('Load sORFDB...')
    sorfdb_df: pl.DataFrame = load_db(args.sorfdb)

    alphabet: pyhmmer.easel.Alphabet = pyhmmer.easel.Alphabet.amino()
    builder: pyhmmer.plan7.Builder = pyhmmer.plan7.Builder(alphabet)
    background: pyhmmer.plan7.Background = pyhmmer.plan7.Background(alphabet)
    hmm: pyhmmer.plan7.HMM

    print('Build HMMs...')
    with (open(args.output.joinpath(f'sorfdb.{SORFDB_VER}.hmm'), 'wb') as hmm_file,
          xopen(args.output.joinpath(f'sorfdb.clusters.{SORFDB_VER}.tsv.gz'), 'w', compresslevel=9, threads=args.threads) as cluster_file):
        cluster_tsv = csv.writer(cluster_file, delimiter='\t', lineterminator='\n')
        cluster_tsv.writerow(('protein', 'cluster'))

        for aln in eio.find_files(args.alignments, r'.*\.afa'):
            with pyhmmer.easel.MSAFile(aln, format='afa') as msa_file:
                msa: pyhmmer.easel.TextMSA = msa_file.read()

                cluster_description: bytes = create_cluster_name(
                    pl.Series('proteins', [protein.decode('UTF-8') for protein in msa.names]),
                    sorfdb_df
                )
                cluster_id: bytes = f"{SORFDB_VER}-{str(aln).split('/')[-1].lstrip('aln.').rstrip('.afa')}".encode('UTF-8')
                msa.name = cluster_id  # cluster ID: id_subgraph.cluster_size.counter

                msa: pyhmmer.easel.DigitalMSA = msa.digitize(alphabet)
                hmm, _, _ = builder.build_msa(msa, background)

                # field names: http://eddylab.org/software/hmmer/Userguide.pdf page 209
                hmm.accession = cluster_id
                if cluster_description is not None and len(cluster_description) > 0:
                    hmm.description = cluster_description

                proteins: list[pyhmmer.easel.DigitalSequence] = msa_proteins_to_digital_seqs(msa.names,
                                                                                             alphabet=alphabet)
                hits: list[namedtuple] = perform_hmmsearch(proteins, hmm, threads=args.threads)

                hmm.cutoffs.gathering = calculate_gathering_cutoff(hits)

                hmm.write(hmm_file)  # , binary=True
                cluster_tsv.writerows(((protein.decode('UTF-8'), cluster_id.decode('UTF-8')) for protein in msa.names))


if __name__ == '__main__':
    main()

# EOF
