#!/usr/bin/env python3

import sys
from pathlib import Path

import extended_io as eio
import constants as const


def high_x_amino_acid_amount(aa: str, x_percentage: float = 0.2) -> bool:
    """
    Check if the protein sequence has a high amount of X amino acids (> 10% or XXX... motif).
    :param aa: protein sequence
    :param x_percentage: sequences with more than x_percentage are to be discarded
    :return: high X amount
    """
    x_count: int = aa.count('X')
    if x_count/len(aa) > x_percentage:
        return True
    return False


def main():
    inverted_fasta: dict[str, str] = dict()

    for file in sys.argv[2:]:
        file_path: Path = Path(str(file).strip()).resolve()
        print('Reading:', str(file_path))
        for header, seq in eio.parse_fasta(file_path):
            if const.MIN_SPROT_LENGTH <= len(seq) <= const.MAX_SPROT_LENGTH:
                if header.startswith('swissprot'):
                    inverted_fasta[seq] = header
                elif not high_x_amino_acid_amount(seq):
                    inverted_fasta[seq] = header

    print('Export combined FAA...')
    faa_outpath: str = str(Path(f"./combined_sorf_db.{sys.argv[1]}.faa.gz"))
    eio.dict_to_fasta({header: aa for aa, header in inverted_fasta.items()}, faa_outpath, threads=1)


if __name__ == '__main__':
    main()

# EOF
