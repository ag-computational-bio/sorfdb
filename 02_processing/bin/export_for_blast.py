#!/usr/bin/env python3

import os
from pathlib import Path
from typing import Any, Generator
from argparse import ArgumentParser

from xopen import xopen

import extended_io as eio
from export_faa import get_faa_entries
import constants as const


def extract_blastx_sprots(json_file: Path) -> Generator[tuple[str, str], Any, None]:
    """
    Add a small protein to the clustering dataset, if its SRV is <= 0.7 (checked with validator.py), is detected by pyrodigal
    and has a canonical start codon. These criteria ensure, that the small protein is highly likely coding.
    :param json_file: blastx genbank json file
    :return: (header, protein)
    """
    for accession_number, bacteria_data in eio.load_json(json_file).items():
        if 'sorfs' not in bacteria_data:
            continue
        for locus_tag, sorf_data in bacteria_data['sorfs'].items():
            if 'uses_sd' in sorf_data and \
                    (sorf_data['orf'].startswith('ATG') or sorf_data['orf'].startswith('GTG') or
                     sorf_data['orf'].startswith('TTG')):
                header: str = f"{accession_number}|{locus_tag}"
                yield header, sorf_data['aa']


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Export all unique small proteins from all databases. '
                                        'Save Entry identifiers to JSON.')
    parser.add_argument('--genbank', '-g', type=Path, help='Input genbank JSON file with sORF data')
    parser.add_argument('--blastx', '-b', type=Path, help='Input BLASTX genbank JSON file with sORF data')
    parser.add_argument('--swissprot', '-s', type=Path, help='Input SwissProtF ASTA file with sORF data')
    parser.add_argument('--uniprot', '-u', type=Path, help='Input UniProt FASTA file with sORF data')
    parser.add_argument('--smprot', '-m', type=Path, help='Input SmProt FASTA file with sORF data')
    parser.add_argument('--output', '-o', type=Path, default=Path('./'), help='Ouput folder path')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    proteins: set[str] = set()
    counter: dict[str, int] = {
        'genbank': 0,
        'blastx': 0,
        'swissprot': 0,
        'uniprot': 0,
        'smprot_combined': 0
    }

    print('Import proteins from genbank...')
    for header, protein in get_faa_entries(args.genbank, unique=False):
        if len(protein) <= const.HMM_MAX_SPROT_LENGTH:
            proteins.add(protein)
            counter['genbank'] += 1

    print('Import proteins from BLASTX genbank...')
    for header, protein in extract_blastx_sprots(args.blastx):
        if len(protein) <= const.HMM_MAX_SPROT_LENGTH:
            proteins.add(protein)
            counter['blastx'] += 1

    print('Import proteins from external DBs...')
    for external_db in (args.swissprot, args.uniprot, args.smprot):
        dbtype: str = str(external_db).split('/')[-1].split('.')[0]
        for header, protein in eio.parse_fasta(external_db):
            if len(protein) <= const.HMM_MAX_SPROT_LENGTH:
                proteins.add(protein)
                counter[dbtype] += 1

    os.makedirs(args.output, exist_ok=True)
    print('Proteins:', len(proteins))
    for dbtype, count in counter.items():
        print(f'\t{dbtype}:', count)

    print('Write FASTA...')
    clustering_faa_path: Path = args.output.joinpath('sprot.clustering.faa.gz')
    with xopen(clustering_faa_path, mode='w', compresslevel=9, threads=args.threads) as f:
        for protein in sorted(list(proteins)):
            f.write(f'>{protein}\n{protein}\n')


if __name__ == '__main__':
    main()

# EOF
