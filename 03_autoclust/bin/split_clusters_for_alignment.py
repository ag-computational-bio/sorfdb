#!/usr/bin/env python3

from pathlib import Path
from argparse import ArgumentParser


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Import a mcl tab file with clusters and export them to fasta files.')
    parser.add_argument('--input', '-i', type=Path, help='Input mcl tab file with clusters.')
    parser.add_argument('--id', '-d', type=str, help='Id of the subgraph.')
    parser.add_argument('--size', '-s', type=int, default=3,
                        help='Minimum sizer of cluster to used for a multiple alignment and HMM.')
    parser.add_argument('--output', '-o', type=Path, default=Path('./'), help='Output folder path')
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    count: int = 0
    ignored: int = 0
    with open(args.input, 'r') as infile:
        for line in infile:
            proteins: list[str] = line.strip().split()
            if len(proteins) >= args.size:
                with open(args.output.joinpath(f'{args.id}.{len(proteins)}.{count}.faa'), 'w') as outfile:
                    for protein in proteins:
                        outfile.write(f'>{protein}\n{protein}\n')
                count += 1
            else:
                ignored += 1

    print(f'Exported clusters: {count}\n'
          f'Ignored clusters:  {ignored}\n'
          f'Total reported:    {count + ignored}')


if __name__ == '__main__':
    main()

# EOF
