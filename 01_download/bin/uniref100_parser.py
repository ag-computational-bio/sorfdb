#!/usr/bin/env python3

import gc
from pathlib import Path
from argparse import ArgumentParser

from xopen import xopen
from lxml import etree as et

import constants as const


def parse_taxon_ids(file_path: Path) -> set[int]:
    """
    Parse the taxon IDs for which eubacteria 2 is the ancestor from a one column TSV file.
    :param file_path: path to taxon id TSV
    :return: bacterial taxon ids
    """
    taxon_ids: set[int] = set()
    with xopen(file_path) as file:
        for i, line in enumerate(file):
            if i > 0:
                taxon_ids.add(int(line.strip().split('\t')[0]))
    return taxon_ids


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Extract sORFs from gbff files.')
    parser.add_argument('--uniref', '-u', type=Path, help='Path to UniRef100 data (XML)')
    parser.add_argument('--taxonomy', '-a', type=str, help='Path to Taxon IDs (TSV).')
    parser.add_argument('--output', '-o', type=Path, default=Path('filtered.fasta.gz'), help='Output path to FASTA file')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    taxon_ids: set[int] = parse_taxon_ids(args.taxonomy)

    with xopen(args.uniref, mode='rb') as fh_xml, xopen(args.output, mode='w', compresslevel=9, threads=args.threads) as out_file:
        for i, (event, elem) in enumerate(et.iterparse(fh_xml, tag='{*}entry')):
            if 'Fragment' not in elem.find('./{*}name').text:  # skip protein fragments
                common_organism = elem.find('./{*}property[@type="common taxon"]')  # source organism
                common_organism = common_organism.get('value') if common_organism is not None else 'Bacteria'

                common_tax_id = elem.find('./{*}property[@type="common taxon ID"]')
                common_tax_id = common_tax_id.get('value') if common_tax_id is not None else 1

                rep_member_dbref = elem.find('./{*}representativeMember/{*}dbReference')
                rep_member_organism = rep_member_dbref.find('./{*}property[@type="source organism"]')  # source organism
                rep_member_organism = rep_member_organism.get('value') if rep_member_organism is not None else 'Bacteria'

                rep_member_tax_id = rep_member_dbref.find('./{*}property[@type="NCBI taxonomy"]')
                rep_member_tax_id = rep_member_tax_id.get('value') if rep_member_tax_id is not None else 1

                if int(common_tax_id) in taxon_ids or int(rep_member_tax_id) in taxon_ids:
                    seq_representative = elem.find('./{*}representativeMember/{*}sequence')
                    seq = seq_representative.text.upper()
                    if len(seq) <= const.MAX_SPROT_LENGTH:
                        uniref100_id = elem.attrib['id']
                        rep_id = rep_member_dbref.attrib['id']
                        name = rep_member_dbref.find('./{*}property[@type="protein name"]')
                        name = name.get('value') if name is not None else 'hypothetical protein'
                        member_count = elem.find('./{*}property[@type="member count"]')
                        member_count = member_count.get('value') if member_count is not None else 1
                        if common_tax_id != 1:
                            tax_id = common_tax_id
                            taxon = common_organism
                        else:
                            tax_id = rep_member_tax_id
                            taxon = rep_member_organism
                        # https://www.uniprot.org/help/fasta-headers
                        # >UniqueIdentifier ClusterName n=Members Tax=TaxonName TaxID=TaxonIdentifier RepID=RepresentativeMember
                        header: str = ' '.join([uniref100_id, name, f'n={member_count}', f'Tax={taxon}',
                                                f'TaxID={tax_id}', f'RepID={rep_id}'])
                        out_file.write(f'>{header}\n{seq}\n')

            elem.clear()
            if i % 100000 == 0:
                del event
                del uniref100_id, rep_id, name, member_count, tax_id, taxon, header
                del seq_representative, seq
                del common_organism, common_tax_id, rep_member_dbref, rep_member_organism, rep_member_tax_id
                gc.collect()


if __name__ == '__main__':
    main()

# EOF
