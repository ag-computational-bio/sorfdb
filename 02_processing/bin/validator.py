#!/usr/bin/env python3

import os
import subprocess as sp
from pathlib import Path
from argparse import ArgumentParser
from typing import Optional, Union

from xopen import xopen

import extended_io as eio
import constants as const


def diamond_search(fasta_path: Path,
                   db_path: Path,
                   out_path: Path,
                   chunks: int = 1,
                   tmpdir: Optional[Path] = None,
                   threads: int = 1) -> dict[str, dict[str, Union[str, float]]]:
    """
    Conduct homology search of hypothetical small proteins against UniProt/smProt db.
    :param fasta_path: Path to FASTA file with hypothetical small proteins
    :param db_path: Path to the combined diamond db
    :param out_path: Output path
    :param chunks: Number of chunks for processing the seed index
    :param tmpdir: Path to diamond tmp directory
    :param threads: Number of threads to use
    :return: Dictionary with hits for small proteins with their DB/homolog name
    """
    # https://github.com/bbuchfink/diamond/discussions/469

    cmd = [
        'diamond',
        'blastp',
        '--db', db_path,
        '--query', fasta_path,
        '--out', out_path,
        '--outfmt', '6', 'qseqid', 'sseqid', 'bitscore', 'evalue',
        '--max-target-seqs', '1',  # single best output
        '--query-cover', '80',  # set mutual coverage to exclude random short small protein to long small protein hits
        '--subject-cover', '80',
        '--matrix', 'BLOSUM62',
        '--ultra-sensitive',
        '--masking', '0',  # -seg
        '--comp-based-stats', '0',  # -composition_base
        '--evalue', '100',
        '--gapped-filter-evalue', '0',
        '--threads', str(threads),
        '--block-size', '0.4',
        '--index-chunks', str(chunks),
        '--compress', '1'
    ]

    if tmpdir:
        # cmd.extend(['--tmpdir', tmpdir])
        pass

    proc = sp.run(
        cmd,
        env=os.environ.copy(),
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if proc.returncode != 0:
        print('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        raise Exception(f'Diamond error! error code: {proc.returncode}')

    hits: dict[str, dict[str, Union[str, float]]] = dict()

    with xopen(out_path, 'r') as fh:
        for line in fh:
            sorf_identifier, cluster_id, bitscore, evalue = line.strip().split('\t')
            hits[sorf_identifier] = {'subject': cluster_id,
                                     'bitscore': float(bitscore)}
    return hits


def extract_valid(hypotheticals: dict[str, dict],
                  srv: dict[str, dict[str, float]],
                  product_references: dict[str, str]) -> dict[str, dict]:
    """
    Return all verified sorf entries from the hypothetical dictionary.
    :param hypotheticals: Dictionary with hypothetical sORFs.
    :param srv: Dictionary containing score ratio values for the sORFs
    :param product_references: Possible product names {{identifiers, ...}: product_name, ...}
    :return: Verified sORFs
    """
    drop_keys: dict[str, list[str]] = dict()

    for accession_number, data in hypotheticals.items():
        if 'sorfs' not in data:
            continue
        for locus_tag, sorf_data in data['sorfs'].items():
            if accession_number in srv and locus_tag in srv[accession_number] and srv[accession_number][locus_tag] >= const.SRV_CUTOFF_VALID:
                # revise hypothetical product names
                identifier: str = f'{accession_number}|{locus_tag}'
                if sorf_data['product'].lower() in const.HYPOTHETICAL_PRODUCT_NAMES:
                    sorf_data['revised_product'] = product_references[identifier]
            else:
                drop_keys.setdefault(accession_number, [])
                drop_keys[accession_number].append(locus_tag)

    for accession_number, locus_tags in drop_keys.items():
        for locus_tag in locus_tags:
            hypotheticals[accession_number]['sorfs'].pop(locus_tag)

    return hypotheticals


def create_diamond_db(faa_path: Path, db_path: Path):
    """
    Create a Diamond DB for the given FAA file.
    :param faa_path: Path to FAA file
    :param db_path: Output path to .dmnd file
    """
    cmd = [
        'diamond',
        'makedb',
        '--in', str(faa_path),
        '--db', str(db_path)
    ]

    proc = sp.run(
        cmd,
        env=os.environ.copy(),
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if proc.returncode != 0:
        print('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        raise Exception(f'Diamond error! error code: {proc.returncode}')


def get_unique_aa(hypotheticals: dict[str, dict]) -> dict[str, str]:
    """
    Create a dictionary for all unique AAs with a header containing all sORFs using this AA.
    :param hypotheticals: Dictionary with hypothetical sORFs
    :return: Dictionary with headers as key and unique AA as value
    """
    aas: dict[str, list[str]] = dict()
    for accession_number, bacteria_data in hypotheticals.items():
        if 'sorfs' in bacteria_data:
            for locus_tag, sorf_data in bacteria_data['sorfs'].items():
                aas.setdefault(sorf_data['aa'], [])
                aas[sorf_data['aa']].append(f"{accession_number}|{locus_tag}")
    return {','.join(headers): aa for aa, headers in aas.items()}


def selfhit_iterator(blastp_data: dict[str, dict[str, Union[str, float]]],
                     sorf_data: dict,
                     out_path: Path,
                     tmpdir: Optional[Path] = None,
                     i: int = 0,
                     threads: int = 1) -> dict[str, dict[str, Union[str, float]]]:
    """
    Return the bitscores for self hits of all hypothetical genes.
    If a hit is no self hit a new Diamond blastp call is executed for the missing hits (recursive).
    :param blastp_data: BLASTP data containing data for all unique AAs self hits or data for missed self hits.
    :param sorf_data: Dictionary with hypothetical sORF data
    :param out_path: Output path
    :param tmpdir: Path to diamond tmp directory
    :param i: Number for the files for missing self hits
    :param threads: Number of threads to use
    :return: data for all self hits
    """
    valid: bool = True
    faa_path: Path = out_path.joinpath(str(i) + '.faa.gz')
    with xopen(faa_path, mode='w', compresslevel=9, threads=threads) as f:
        for headers, data in blastp_data.items():
            if headers != data['subject']:
                valid = False
                first_header: str = headers.split(',')[0]  # Needed to get original aa
                accession_number: str = first_header.split('|')[0]
                locus_tag: str = '|'.join(first_header.split('|')[1:])
                f.write(f">{headers}\n{sorf_data[accession_number]['sorfs'][locus_tag]['aa']}\n")

    if valid:
        return blastp_data

    db_path: Path = out_path.joinpath(str(i) + '.dmnd')
    create_diamond_db(faa_path, db_path)
    db_out_path: Path = out_path.joinpath(str(i) + '.tsv.gz')
    i += 1

    blastp_data.update(selfhit_iterator(diamond_search(faa_path,
                                                       db_path,
                                                       db_out_path,
                                                       tmpdir=tmpdir,
                                                       threads=threads),
                                        sorf_data=sorf_data,
                                        out_path=out_path,
                                        tmpdir=tmpdir,
                                        i=i,
                                        threads=threads))
    return blastp_data


def parse_self_hits(self_hits: dict[str, dict[str, Union[str, float]]]) -> dict[str, dict[str, float]]:
    """
    Return the bitscore for each sorf in a dictionary accassible via {accession_number: {locus_tag: bitscore}}.
    :param self_hits: BLASTP data for all self hits
    :return: bitscores mapped to their accession_number and locus_tag
    """
    self_hit_counter: int = 0
    duplicate_counter: int = 0
    self_bitscores: dict[str, dict[str, float]] = dict()
    for headers, data in self_hits.items():
        header_ls: list[str] = headers.split(',')
        for header in header_ls:
            self_hit_counter += 1
            tmp: list[str] = header.split('|')
            accession_number: str = tmp[0]
            locus_tag: str = '|'.join(tmp[1:])
            self_bitscores.setdefault(accession_number, dict())
            # if locus_tag in self_bitscores[accession_number]:  # DEBUG
            #     print(accession_number, locus_tag)
            #     duplicate_counter += 1
            self_bitscores[accession_number][locus_tag] = data['bitscore']

    assert self_hit_counter == sum([len(x) for x in self_bitscores.values()]), \
        f"Number of self-hits={self_hit_counter} != " \
        f"Number of bitscores={sum([len(x) for x in self_bitscores.values()])}\n" \
        f"DUPLICATES: {duplicate_counter}"

    return self_bitscores


def calculate_srv(self_bitscores: dict[str, dict[str, float]],
                  homolog_hits: dict[str, dict[str, Union[str, float]]],
                  data: dict) -> dict[str, dict[str, float]]:
    """
    Calculate the Score Ratio Values for all sORFs with self-hits and hits in the DB.
    :param self_bitscores: Dictionary with self-hit bitscores {accession_number: {locus_tag: bitscore}}
    :param homolog_hits: BLASTP hit data of the hypothetical sORFs
    :param data: hypothetical sORF data
    :return: Dictionary with SRVs {accession_number: {locus_tag: SRV}}
    """
    print('\tCalculate Score Ratio Values...')
    srv: dict[str, dict[str, float]] = dict()

    found: int = 0
    missed: int = 0
    for accession_number, data in data.items():
        if accession_number == 'metadata':
            continue
        srv.setdefault(accession_number, dict())
        for locus_tag, sorf_data in data['sorfs'].items():
            identifier: str = f"{accession_number}|{locus_tag}"
            if (accession_number in self_bitscores and locus_tag in self_bitscores[accession_number]) and identifier in homolog_hits:
                # SRV = bitscore(a->b)/bitscore(a->a)
                self_hit: float = self_bitscores[accession_number][locus_tag]
                homolog_hit: float = homolog_hits[identifier]['bitscore']
                srv[accession_number][locus_tag] = round(homolog_hit/self_hit, 4)
                found += 1
            else:
                missed += 1

    print(f"\t\tFound: {found}\n\t\tMissed: {missed}")
    return srv


def get_possible_product_name(header: str) -> str:
    """
    Get possible product name from UniProt/SwissProt subject id.
    :param header: Name of the best hit cluster (fasta header)
    :return: Possible product name
    """
    # https://www.uniprot.org/help/fasta-headers
    # >db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
    # >UniqueIdentifier ClusterName n=Members Tax=TaxonName TaxID=TaxonIdentifier RepID=RepresentativeMember
    header_as_list: list[str] = header.split()
    stop: int = 1
    if header_as_list[0].startswith('swissprot|UniRef') or header_as_list[0].startswith('uniprot|UniRef'):
        ref_id: str = header_as_list[0].split('|')[1]

        for i, x in enumerate(header_as_list):
            if x.startswith('n='):
                stop = i
                break
        product: str = ' '.join(header_as_list[1:stop] + [f'({ref_id})'])
    else:
        try:
            unique_identifier, entry_name = header_as_list[0].split('|')[-2:]
        except ValueError:
            print(print(f'Could not split header to extract product:\n{header}'))
            exit(1)

        ref_id: str = f'{unique_identifier}|{entry_name}'
        for i, x in enumerate(header_as_list):
            if x.startswith('OS='):
                stop = i
                break
        product: str = ' '.join(header_as_list[1:stop] + [f'({ref_id})'])
    return product


def load_product_names(database_path: Path, hits: dict[str, dict[str, Union[str, float]]]) -> dict[str, str]:
    """
    Load the database FASTA file and extract the product names from the headers.
    :param database_path: Path to the FASTA file from which the diamond db was created
    :param hits: Best diamond blastp hits against the known small protein database
    :return: Possible product names {identifier: product_name, ...}
    """
    protein_products: dict[str, str] = dict()
    for faa in ('swissprot.faa.gz', 'uniprot.filtered.faa.gz'):
        with xopen(database_path.joinpath(faa).resolve(), mode='r') as f:
            for line in f:
                if line.startswith('>'):
                    line = line.strip().lstrip('>')
                    protein_products[line.split()[0]] = get_possible_product_name(line)  # id: product

    product_references: dict[str, str] = dict()
    for identifiers, blast_hit in hits.items():
        if blast_hit['subject'].startswith('swissprot') or blast_hit['subject'].startswith('uniprot'):
            product = protein_products[blast_hit['subject']]
        else:
            product = blast_hit['subject']  # smProt
        for identifier in identifiers.split(','):
            product_references[identifier] = product

    return product_references


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Validate hypothetical sORFs with diamond blastp against UniProt and smProt.')
    parser.add_argument('--hypotheticals', type=Path, help='Path to JSON file with hypothetical sORFs')
    parser.add_argument('--fasta', type=Path, help='Path to FASTA file with hypothetical sORFs')
    parser.add_argument('--db', type=Path, help='Path to the protein database.')
    parser.add_argument('--tmpdir', type=Path, default=None, help='Diamond tmpdir')
    parser.add_argument('--source', type=str, default='genbank', help='Source (default=genbank)')
    parser.add_argument('--output', type=Path, default=Path('./'), help='Output path')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads')
    return parser.parse_args()


def main():
    args = parse_arguments()
    hypotheticals_path, fasta_path, source = args.hypotheticals, args.fasta, args.source
    db_path: Path = args.db
    tmpdir: Optional[Path] = args.tmpdir
    out_path: Path = args.output
    threads: int = args.threads

    print('Homolog search...')
    diamond_output_path: Path = out_path.joinpath('diamond.sorfs.tsv.gz')
    best_homolog_hits: dict[str, dict[str, Union[str, float]]] = diamond_search(fasta_path,
                                                                                db_path,
                                                                                diamond_output_path,
                                                                                tmpdir=tmpdir,
                                                                                threads=threads)

    print('Load data...')
    hypothetical_sorf_data: dict = eio.load_json(hypotheticals_path)
    count_hypotheticals: int = sum([len(x['sorfs']) for x in hypothetical_sorf_data.values() if 'sorfs' in x])

    print('Determine Score Ratio Values...')
    print('\tWrite unique AAs...')
    unique_hypothetical_faa: Path = out_path.joinpath('unique.hypotheticals.faa.gz')
    eio.dict_to_fasta(get_unique_aa(hypothetical_sorf_data),
                      unique_hypothetical_faa,
                      threads=threads)
    print('\tCreate self-hit DB...')
    unique_hypothetical_db: Path = out_path.joinpath('unique.hypotheticals.dmnd')
    create_diamond_db(unique_hypothetical_faa,
                      unique_hypothetical_db)
    print('\tGet self-hits...')
    unique_hypothetical_output_path: Path = out_path.joinpath('diamond.selfhit.tsv.gz')
    srv: dict[str, dict[str, float]]
    srv = calculate_srv(parse_self_hits(selfhit_iterator(diamond_search(unique_hypothetical_faa,
                                                                        unique_hypothetical_db,
                                                                        unique_hypothetical_output_path,
                                                                        tmpdir=tmpdir,
                                                                        threads=threads),
                                                         sorf_data=hypothetical_sorf_data,
                                                         out_path=out_path,
                                                         tmpdir=tmpdir,
                                                         threads=threads)),
                        best_homolog_hits,
                        hypothetical_sorf_data)

    count_valid_srv: int = len([1 for x in srv.values() for y in x.values() if y >= const.SRV_CUTOFF_VALID])

    # Dump Score Ratio Values for debugging
    eio.export_to_json(srv, out_path.joinpath('srv.json.gz'), threads=threads)

    db_base_path: Path = Path('/'.join(str(args.db).split('/')[:-1]))  # get db base path
    product_references: dict[str, str] = load_product_names(db_base_path,
                                                            hits=best_homolog_hits)
    del best_homolog_hits
    print('Extract verified sORFs...')
    hypothetical_sorf_data = extract_valid(hypotheticals=hypothetical_sorf_data,
                                           srv=srv,
                                           product_references=product_references)

    count_verified: int = sum([len(x['sorfs']) for x in hypothetical_sorf_data.values() if 'sorfs' in x])
    assert count_verified == count_valid_srv, f"Number of verified sORFs={count_verified} does not match expected number={count_valid_srv}"
    print('\tHypothetical sORFs:', count_hypotheticals)
    print('\tVerified sORFs:', count_verified)

    print('Export verified sORFs to JSON...')
    result_output_path: Path = out_path.joinpath(f"verified.{source}.json")

    eio.export_to_json(hypothetical_sorf_data, result_output_path, threads=threads)


if __name__ == '__main__':
    main()

# EOF
