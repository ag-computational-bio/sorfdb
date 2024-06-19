#!/usr/bin/env python3

from pathlib import Path
from argparse import ArgumentParser

from xopen import xopen

import extended_io as eio


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Extract sORFs from gbff files.')
    parser.add_argument('--input', '-i', type=Path, help='Path to UniProt/SwissProt data (FASTA)')
    parser.add_argument('--uniref', '-u', type=Path, help='Path to UniRef100 data (FASTA)')
    parser.add_argument('--prefix', '-p', type=str, help='Header prefix')
    parser.add_argument('--output', '-o', type=Path, default=Path('./out.fasta.gz'), help='Output path to FASTA file')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    query_data: dict[str, str] = {protein: header for header, protein in eio.parse_fasta(args.input)}
    with xopen(args.output, 'w', compresslevel=9, threads=args.threads) as outfile:
        for header, protein in eio.parse_fasta(args.uniref):
            if protein in query_data:
                outfile.write(f'>{args.prefix}|{header}\n{protein}\n')
                query_data.pop(protein)
        print('Missed:', len(query_data), 'entries in', args.prefix)
        for protein, header in query_data.items():
            outfile.write(f'>{args.prefix}|{header}\n{protein}\n')


if __name__ == '__main__':
    main()

# EOF
