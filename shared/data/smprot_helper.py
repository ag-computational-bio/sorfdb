#!/usr/bin/env python3

import sys
from pathlib import Path

import extended_io as eio


def parse_smprot(file_path: str, db: str) -> dict[str, str]:
    """
    Parse the smProt database files (txt/tsv format.)
    :param file_path: Path to smProt database file [txt|tsv]
    :param db: smProt sub-db
    :return: Dictionary with smProt small proteins.
    """
    smprot_dict: dict[str, str] = dict()

    with open(file_path, 'r') as f:
        for line in f.readlines():
            # SmProt_ID RNASeq ProteinSequence ProteinLength Chr Start Stop Strand Blocks StartCodon PhyloCSF_Mean
            # minRiboP minTISP
            split_line: list[str] = line.rstrip('\n').split()
            if not split_line[0].startswith('#') and split_line[2] != 'NA':
                genome: str = split_line[4]
                smprot_id: str = split_line[0]
                aa: str = split_line[2].rstrip('*')
                smprot_dict[f"{db}|{genome}|{smprot_id}"] = aa
    return smprot_dict


def main():
    _, smprot_literature, smprot_kdb, smprot_ribo_seq = sys.argv
    out_path = Path('./')

    for db_name, data in zip(('smprot_literature', 'smprot_kdb', 'smprot_ribo_seq'),
                             (smprot_literature, smprot_kdb, smprot_ribo_seq)):
        faa_path: str = str(out_path.joinpath(f"{db_name}.faa.gz"))
        eio.dict_to_fasta(parse_smprot(data, db_name), faa_path)


if __name__ == '__main__':
    main()

# EOF
