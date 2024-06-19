#!/usr/bin/env python3

from pathlib import Path
from argparse import ArgumentParser

import ijson
from xopen import xopen


def get_faa_entries(json_file: Path, unique: bool = False):
    unique_aa: set[str] = set()
    with xopen(json_file, mode='rb') as f:
        parser = ijson.parse(f, use_float=True)
        accession_number, locus_tag, aa = '', '', ''
        for prefix, event, value in parser:
            if len(prefix) == 0 and event == 'map_key':
                accession_number: str = value
            elif event == 'map_key' and prefix == f"{accession_number}.sorfs":
                locus_tag: str = value
            elif prefix == f"{accession_number}.sorfs.{locus_tag}.aa":
                aa: str = value
                if unique:
                    if aa not in unique_aa:
                        header: str = f"{accession_number}|{locus_tag}"
                        yield header, aa
                        unique_aa.add(aa)
                else:
                    header: str = f"{accession_number}|{locus_tag}"
                    yield header, aa
                locus_tag, aa = '', ''


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Output the upstream sequences of the given json data set.')
    parser.add_argument('--input', '-i', type=Path, help='Input JSON file with sORF data')
    parser.add_argument('--output', '-o', type=Path, default=Path('./'), help='Output file path')
    parser.add_argument('--unique', '-u', action='store_true', help='Export only unique sequences.')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    return args.input, args.output, args.unique, args.threads


def main():
    in_file, faa_path, unique, threads = parse_arguments()
    print(f'Get unique={unique} faa sequences.')

    ijson.get_backend('yajl2_cffi')

    with xopen(faa_path, mode='w', compresslevel=9, threads=threads) as f:
        for header, aa in get_faa_entries(json_file=in_file, unique=unique):
            f.write(f">{header}\n{aa}\n")


if __name__ == '__main__':
    main()

# EOF
