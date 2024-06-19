#!/usr/bin/env python3

import os
import csv
import json
import tarfile
from io import BytesIO
from time import time
from pathlib import Path
from argparse import ArgumentParser
from typing import Generator, Optional, Union

import pyhmmer
import polars as pl
from xopen import xopen
from peptides import Peptide
from Bio.SeqUtils.ProtParam import ProteinAnalysis

import constants as const
import extended_io as eio
from hmm import perform_hmmsearch, reformat_results


new_phyla: dict[str, str] = {
    'Acidobacteria': 'Acidobacteriota',
    'Actinobacteria': 'Actinomycetota',
    'Aquificae': 'Aquificota',
    'Armatimonadetes': 'Armatimonadota',
    'Atribacter': 'Atribacterota',
    'Firmicutes': 'Bacillota',
    'Bacteroidetes': 'Bacteroidota',
    'Balneolaeota': 'Balneolota',
    'Bdellovibrio': 'Bdellovibrionota',
    'Caldiserica': 'Caldisericota',
    'Calditrichaeota': 'Calditrichota',
    'Campylobacteria': 'Campylobacterota',
    'Epsilonbacteraeota': 'Campylobacterota',
    'Chlamydiae': 'Chlamydiota',
    'Chlorobi': 'Chlorobiota',
    'Chloroflexi': 'Chloroflexota',
    'Chrysiogenetes': 'Chrysiogenota',
    'Coprothermobacter': 'Coprothermobacterota',
    'Deferribacteres': 'Deferribacterota',
    'Deinococcus': 'Deinococcota',
    'Deinococcus-Thermus': 'Deinococcota',
    'Dictyoglomi': 'Dictyoglomota',
    'Elusimicrobia': 'Elusimicrobiota',
    'Fibrobacteres': 'Fibrobacterota',
    'Fusobacteria': 'Fusobacteriota',
    'Gemmatimonadetes': 'Gemmatimonadota',
    'Ignavibacteriae': 'Ignavibacteriota',
    'Kiritimatiellaeota': 'Kiritimatiellota',
    'Lentisphaerae': 'Lentisphaerota',
    'Tenericutes': 'Mycoplasmatota',
    'Myxococcus': 'Myxococcota',
    'Thaumarchaeota': 'Nitrososphaerota',
    'Nitrospinae': 'Nitrospinota',
    'Nitrospirae': 'Nitrospirota',
    'Planctomycetes': 'Planctomycetota',
    'Proteobacteria': 'Pseudomonadota',
    'Rhodothermaeota': 'Rhodothermota',
    'Spirochaetes': 'Spirochaetota',
    'Synergistetes': 'Synergistota',
    'Thermodesulfobacteria': 'Thermodesulfobacteriota',
    'Desulfobacterota': 'Thermodesulfobacteriota',
    'Thermomicrobia': 'Thermomicrobiota',
    'Crenarchaeota': 'Thermoproteota',
    'Thermotogae': 'Thermotogota',
    'Verrucomicrobia': 'Verrucomicrobiota'
}


def parse_taxon_ids(file_path: Path) -> dict[int, dict[str, str]]:
    """
    Parse the taxon IDs for which eubacteria 2 is the ancestor from a TSV file.
    Columns:
    TaxID Scientific_name Lineage
    :param file_path: path to taxon id TSV
    :return: bacterial taxa by taxid
    """
    global new_phyla
    taxon_keys: tuple[str, str, str, str, str, str] = ('phylum', 'class', 'order', 'family', 'genus', 'species')
    taxon_by_taxid: dict[int, dict[str, str]] = {
        2: {'phylum': 'Bacteria', 'class': '', 'order': '', 'family': '', 'genus': '', 'species': '', 'strain': ''}
    }
    taxid: int
    scientific_name: str
    lineage_str: str
    with xopen(file_path) as file:
        _ = file.readline().split('\t')
        for line in file:
            taxid, strain, lineage_str = line.strip().split('\t')
            taxid = int(taxid)
            lineage: list[str] = [elem.strip() for elem in lineage_str.split(', ')]
            try:
                index: int = last_index(lineage, 'Bacteria')
                lineage = lineage[index+1:]
            except ValueError:
                lineage = []
            try:
                if lineage[0].endswith('group'):
                    lineage = lineage[1:]
            except IndexError:
                pass  # no taxonomy available

            taxon_by_taxid[taxid] = {'strain': strain}
            for i, taxon_key in enumerate(taxon_keys):
                try:
                    if i == 0 and lineage[i] in new_phyla:
                        taxon_by_taxid[taxid][taxon_key] = new_phyla[lineage[i]]
                    else:
                        taxon_by_taxid[taxid][taxon_key] = lineage[i]
                except IndexError:
                    taxon_by_taxid[taxid][taxon_key] = ''
    return taxon_by_taxid


def last_index(listlike: list, search_element) -> int:
    """
    Return the index of the last occurrence of an element in a list (or another iterable).
    If elemnet not in list rise a ValueError.
    :param listlike: list
    :param search_element: element to look up
    :return: last index of the search element
    """
    index: int = -1
    for i, element in enumerate(listlike):
        if element == search_element:
            index = i
    if index != -1:
        return index
    else:
        raise ValueError(f'Element {search_element} not in list.')


def revise_taxonomy(taxonomy: tuple[str, list[str]]) -> dict[str, str]:
    """
    Standardize the taxonomy of genbank entries. Handle missing fields. Handle wrong entries.
    :param taxonomy: Taxonomy (species/strain (phylum, ...))
    :return: Dictionary with taxonmy {phylum: ..., ...}
    """
    taxonomy_revised: dict[str, str] = dict()
    tree: list[str] = taxonomy[1]

    for i, tax_key in enumerate(('phylum', 'class', 'order', 'family', 'genus', 'species', 'strain')):
        try:
            taxonomy_revised[tax_key] = tree[i + 1]  # +1 ro skip kingdom entry
        except IndexError:  # set empty taxonomy entries
            taxonomy_revised[tax_key] = ''

    if taxonomy[0].count(' ') > 1:  # if strain information are available put them in a separate field
        taxonomy_revised['species'] = ' '.join(taxonomy[0].split()[0:2])
        taxonomy_revised['strain'] = ' '.join(taxonomy[0].split()[2:])
    else:
        taxonomy_revised['species'] = taxonomy[0]
        taxonomy_revised['strain'] = ''
    return taxonomy_revised


def advanced_gravy(prot_analysis: ProteinAnalysis, protein: str) -> Optional[float]:
    """
    Handle Selenocystein and Pyrolysin for protein gravy calculation. Replace U -> C and O -> K.
    :param prot_analysis: ProteinAnalysis object
    :param protein: protein sequence
    :return: gravy
    """
    try:
        return round(prot_analysis.gravy(), 4)
    except KeyError as e:
        if e == 'U' or e == 'O':
            prot_analysis: ProteinAnalysis = ProteinAnalysis(protein.replace('U', 'C').replace('O', 'K'))
            return round(prot_analysis.gravy(), 4)
        else:
            return None


def advanced_instability(prot_analysis: ProteinAnalysis, protein: str) -> Optional[float]:
    """
    Handle Selenocystein and Pyrolysin for protein instability_index calculation. Replace U -> C and O -> K.
    :param prot_analysis: ProteinAnalysis object
    :param protein: protein sequence
    :return: instability_index
    """
    try:
        return round(prot_analysis.instability_index(), 4)
    except KeyError as e:
        if e == 'U' or e == 'O':
            prot_analysis: ProteinAnalysis = ProteinAnalysis(protein.replace('U', 'C').replace('O', 'K'))
            return round(prot_analysis.instability_index(), 4)
        else:
            return None


def advanced_molecular_weight(prot_analysis: ProteinAnalysis, protein: str) -> Optional[float]:
    """
    Handle Selenocystein and Pyrolysin for protein molecular_weight calculation. Replace U -> C and O -> K.
    :param prot_analysis: ProteinAnalysis object
    :param protein: protein sequence
    :return: molecular_weight
    """
    try:
        return round(prot_analysis.molecular_weight(), 4)
    except ValueError as e:
        if e == 'U' or e == 'O':
            prot_analysis: ProteinAnalysis = ProteinAnalysis(protein.replace('U', 'C').replace('O', 'K'))
            return round(prot_analysis.molecular_weight(), 4)
        else:
            return None


def genbank_extract_db_fields(genbank_data: dict) -> Generator[dict[str, Union[str, float]], None, None]:

    global new_phyla

    for accession_number, data in genbank_data.items():
        if 'sorfs' not in data:
            continue
        taxonomy: dict[str, str] = revise_taxonomy(data['taxonomy'])
        for protein_id, sorf_data in data['sorfs'].items():
            prot_analysis: ProteinAnalysis = ProteinAnalysis(sorf_data['aa'])
            peptide = Peptide(sorf_data['aa'])

            datapoint: dict[str, Union[str, float]] = {
                'id': f'GenBank|{accession_number}|{protein_id}',
                'source': 'GenBank',
                'assembly': data['assembly'],
                'accession': accession_number,
                'protein-id': protein_id,
                'uid': '',
                'entry-name': '',
                'phylum': new_phyla[taxonomy['phylum']] if taxonomy['phylum'] in new_phyla else taxonomy['phylum'],
                'class': taxonomy['class'],
                'order': taxonomy['order'],
                'family': taxonomy['family'],
                'genus': taxonomy['genus'],
                'species': taxonomy['species'],
                'strain': taxonomy['strain'],
                'sorf': sorf_data['orf'],
                'slen': len(sorf_data['orf']),
                'start-codon': sorf_data['orf'][:3],
                'protein': sorf_data['aa'],
                'plen': len(sorf_data['aa']),
                'product': sorf_data['revised_product'] if 'revised_product' in sorf_data else sorf_data['product'],
                'rbs': sorf_data.get('uses_sd', 0),
                'pfam-hits': [{'name': hit, 'evalue': sorf_data['pfam-evalue'][i]} for i, hit in
                              enumerate(sorf_data['pfam-subject'])] if 'pfam-subject' in sorf_data else [],
                'gravy': advanced_gravy(prot_analysis, sorf_data['aa']),  # hydrophobicity
                'aromaticity': round(prot_analysis.aromaticity(), 4),
                'molecular-weight': advanced_molecular_weight(prot_analysis, sorf_data['aa']),
                'instability': advanced_instability(prot_analysis, sorf_data['aa']),
                'isoelectric-point': round(prot_analysis.isoelectric_point(), 4),
                'aliphatic-index': round(peptide.aliphatic_index(), 4),
                'boman': round(peptide.boman(), 4)
            }
            yield datapoint


def extract_blastx_genbank(json_file: Path) -> dict:
    """
    Add a small protein to the sORFDB dataset, if its SRV is <= 0.7 (checked with validator.py), is detected by pyrodigal
    and has a canonical start codon. These criteria ensure, that the small protein is highly likely coding.
    :param json_file: blastx genbank json file
    :return: filtered blastx genbank json file
    """
    high_quality_bastx: dict = dict()
    for accession_number, bacteria_data in eio.load_json(json_file).items():
        if 'sorfs' not in bacteria_data:
            continue
        for locus_tag, sorf_data in bacteria_data['sorfs'].items():
            if 'uses_sd' in sorf_data and \
                    (sorf_data['orf'].startswith('ATG') or sorf_data['orf'].startswith('GTG') or
                     sorf_data['orf'].startswith('TTG')):
                high_quality_bastx.setdefault(accession_number,
                                              {'assembly': bacteria_data['assembly'],
                                               'taxonomy': bacteria_data['taxonomy'],
                                               'sorfs': dict()}
                                              )
                high_quality_bastx[accession_number]['sorfs'][locus_tag] = sorf_data
    return high_quality_bastx


def get_possible_product_name(fasta_header: str, source: str) -> tuple[str, str, str, int]:
    """
    Get possible product name from UniProt/SwissProt/SmProt subject id.
    :param fasta_header: Name of the best hit cluster
    :param source: source database
    :return: Possible product name
    """
    taxid: int = 2

    # https://www.uniprot.org/help/fasta-headers
    header_as_list: list[str] = fasta_header.split()
    try:
        unique_identifier, entry_name = header_as_list[0].split('|')[1:]  # SwissProt and UniProt proteinIDs
    except ValueError:
        print(f'Splitting error for:\n{fasta_header}')
        exit(1)

    stop: int = 1
    for i, x in enumerate(header_as_list):
        if x.startswith('OS='):
            stop = i
        elif x.startswith('OX='):
            taxid = int(x.lstrip('OX='))
            break
    else:
        print(f'Could not get product and TaxID for:\n{fasta_header}')
        exit(1)
    product: str = ' '.join(header_as_list[1:stop])
    return product, unique_identifier, entry_name, taxid


def fasta_extract_db_fields(fasta_path: Path,
                            pfam_path: Path,
                            taxonomy_by_id:  dict[int, dict[str, str]],
                            threads: int = 1) -> Generator[dict[str, Union[str, float]], None, None]:
    """
    Create database datapoints (JSON) for external DBs (SwissProt/UniProt/SmProt).
    :param fasta_path: path to FASTA file
    :param pfam_path: path to pfam HMM file
    :param taxonomy_by_id: dictionary with taxid keys and taxonomy as value
    :param threads: number of threads
    :return:
    """
    source_db: str = str(fasta_path).split("/")[-1].split(".")[0]
    if source_db == 'swissprot':
        source: str = 'SwissProt'
    elif source_db == 'uniprot':
        source: str = 'UniProt'
    else:
        raise ValueError(f'Unknown database={source_db}')

    print(f'Export {source} sORFs...')

    pfam = reformat_results(perform_hmmsearch(fasta_proteins_to_digital_seqs(fasta_path),
                                              str(pfam_path),
                                              cutoff='gathering',
                                              best_hit=False,
                                              threads=threads),
                            evalue=1E-5,
                            prefix='pfam')

    for header, seq in eio.parse_fasta(fasta_path):
        product, unique_identifier, entry_name, taxid = get_possible_product_name(header, source)
        try:
            taxonomy: dict[str, str] = taxonomy_by_id[taxid]
        except KeyError:
            print(f'TaxID error for:\n{header}')
            exit(1)

        prot_analysis: ProteinAnalysis = ProteinAnalysis(seq)
        peptide = Peptide(seq)
        pfam_entries: dict[str, list[float | str] | str] = get_pfam_entries(header, pfam)
        pfam_hits = []
        if len(pfam_entries['pfam-subject']) > 0:
            try:
                pfam_hits = [{'name': hit, 'evalue': pfam_entries['pfam-evalue'][i]} for i,  hit in enumerate(pfam_entries['pfam-subject'])]
            except KeyError:
                pass

        datapoint: dict[str, Union[str, float]] = {
            'id': f'{source}|{unique_identifier}|{entry_name}',
            'source': source,
            'assembly': '',
            'accession': '',
            'protein-id': '',
            'uid': unique_identifier,
            'entry-name': entry_name,
            'phylum': taxonomy['phylum'],
            'class': taxonomy['class'],
            'order': taxonomy['order'],
            'family': taxonomy['family'],
            'genus': taxonomy['genus'],
            'species': taxonomy['species'],
            'strain': taxonomy['strain'],
            'sorf': '',
            'slen': len(seq) * 3 + 3,
            'start-codon': '',
            'protein': seq,
            'plen': len(seq),
            'product': product,
            'rbs': None,
            'pfam-hits': pfam_hits,
            'gravy': advanced_gravy(prot_analysis, seq),  # hydrophobicity
            'aromaticity': round(prot_analysis.aromaticity(), 4),
            'molecular-weight': advanced_molecular_weight(prot_analysis, seq),
            'instability': advanced_instability(prot_analysis, seq),
            'isoelectric-point': round(prot_analysis.isoelectric_point(), 4),
            'aliphatic-index': round(peptide.aliphatic_index(), 4),
            'boman': round(peptide.boman(), 4)
        }
        yield datapoint


def smprot_extract_db_fields(db_path: Path,
                             pfam_path: Path,
                             threads: int = 1) -> Generator[dict[str, Union[str, float]], None, None]:
    """
    Create database datapoints (JSON) for external DBs (SwissProt/UniProt/SmProt).
    :param db_path: path to SmProt combined TSV file
    :param pfam_path: path to pfam HMM file
    :param threads: number of threads
    :return:
    """
    df = pl.scan_csv(
        db_path,
        separator='\t'
    ).select(
        pl.col('Chr'),
        pl.col('#SmProt_ID'),
        pl.col('RNASeq'),
        pl.col('ProteinSequence')
    ).rename(
        {
            'Chr': 'accession',
            '#SmProt_ID': 'id',
            'RNASeq': 'sorf',
            'ProteinSequence': 'protein'
        }
    ).with_columns(
        pl.col('protein').str.rstrip('*'),
        pl.col('sorf').str.replace('NA', '', literal=True)
    ).unique().collect()

    sequences: list[pyhmmer.easel.DigitalSequence] = []
    for accession, protid, aa in df.select(pl.col('accession'), pl.col('id'), pl.col('protein')).iter_rows():
        sequences.append(
            pyhmmer.easel.TextSequence(sequence=aa,
                                       accession=bytes(accession, 'utf-8'),
                                       name=bytes(protid, 'utf-8')).digitize(
                pyhmmer.easel.Alphabet.amino())
        )

    pfam = reformat_results(perform_hmmsearch(sequences,
                                              str(pfam_path),
                                              cutoff='gathering',
                                              best_hit=False,
                                              threads=threads),
                            evalue=1E-5,
                            prefix='pfam')
    del sequences

    for accession, protid, sorf, protseq in df.iter_rows():
        prot_analysis: ProteinAnalysis = ProteinAnalysis(protseq)
        peptide = Peptide(protseq)
        try:
            pfam_hits = [{'name': hit, 'evalue': pfam[accession][protid]['pfam-evalue'][i]} for i,  hit in enumerate(pfam[accession][protid]['pfam-subject'])]
        except KeyError:
            pfam_hits = []

        datapoint: dict[str, Union[str, float]] = {
            'id': f'SmProt|{accession}|{protid}',
            'source': 'SmProt',
            'assembly': '',
            'accession': accession,
            'protein-id': protid,
            'uid': '',
            'entry-name': '',
            'phylum': 'Pseudomonadota',
            'class': 'Gammaproteobacteria',
            'order': 'Enterobacterales',
            'family': 'Enterobacteriaceae',
            'genus': 'Escherichia',
            'species': 'Escherichia coli',
            'strain': '',
            'sorf': sorf,
            'slen': len(sorf) if len(sorf) > 0 else len(protseq) * 3 + 3,
            'start-codon': sorf[:3],
            'protein': protseq,
            'plen': len(protseq),
            'product': 'hypothetical protein',
            'rbs': None,
            'pfam-hits': pfam_hits,
            'gravy': advanced_gravy(prot_analysis, protseq),  # hydrophobicity
            'aromaticity': round(prot_analysis.aromaticity(), 4),
            'molecular-weight': advanced_molecular_weight(prot_analysis, protseq),
            'instability': advanced_instability(prot_analysis, protseq),
            'isoelectric-point': round(prot_analysis.isoelectric_point(), 4),
            'aliphatic-index': round(peptide.aliphatic_index(), 4),
            'boman': round(peptide.boman(), 4)
        }
        yield datapoint


def fasta_proteins_to_digital_seqs(fasta_file: Path) -> list[pyhmmer.easel.DigitalSequence]:
    """
    Convert protein sequences of a FASTA file to a list of digital pyhmmer Sequences.
    :param fasta_file: path to FASTA file
    :return: list of pyhmmer DigitalSequences from the sORF proteins
    """
    sequences: list[pyhmmer.easel.DigitalSequence] = []
    for header, aa in eio.parse_fasta(fasta_file):
        header = header.strip().split()[0]
        header_ls = header.split('|')
        sequences.append(
            pyhmmer.easel.TextSequence(sequence=aa,
                                       accession=bytes(header_ls[0], 'utf-8'),
                                       name=bytes('|'.join(header_ls[1:]), 'utf-8')).digitize(
                pyhmmer.easel.Alphabet.amino())
        )
    return sequences


def get_pfam_entries(identifier: str,
                     pfam_hits: dict[str, dict[str, dict[str, list[Union[float, str]]]]]) -> dict[str, Union[list[Union[float, str]], str]]:
    """
    Return the Pfam hits for an external db protein from the pfam_hits.
    :param identifier: header of the protein (accession, protein_id)
    :param pfam_hits: pfam hits
    :return: pfam hits for the given protein
    """
    header_ls: list[str] = identifier.strip().split()[0].split('|')
    accession: str = header_ls[0]
    name: str = '|'.join(header_ls[1:])
    try:
        return pfam_hits[accession][name]
    except KeyError:
        return {
            'pfam-subject': '',
            'pfam-evalue': ''
        }


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Export all sORFs to JSON datapoints for Elasticsearch and to TSV.')
    parser.add_argument('--genbank', '-g', type=Path, help='Input genbank JSON file with sORF data')
    parser.add_argument('--blastx', '-b', type=Path, help='Input blastx genbank JSON file with sORF data')
    parser.add_argument('--swissprot', '-s', type=Path, help='Input SwissProtF ASTA file with sORF data')
    parser.add_argument('--uniprot', '-u', type=Path, help='Input UniProt FASTA file with sORF data')
    parser.add_argument('--smprot', '-m', type=Path, help='Input SmProt FASTA file with sORF data')
    parser.add_argument('--taxonomy', '-x', type=Path, help='Input bacterial taxonomy TSV with TaxID, '
                                                            'scientific name and lineage')
    parser.add_argument('--pfam', '-p', type=Path, help='Input Pfam hmm')
    parser.add_argument('--output', '-o', type=Path, default=Path('./'), help='Output folder path')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    os.makedirs(args.output, exist_ok=True)

    db_fields: dict[str, list] = {
        'id': [],
        'source': [],
        'assembly': [],
        'accession': [],
        'protein-id': [],
        'uid': [],
        'entry-name': [],
        'phylum': [],
        'class': [],
        'order': [],
        'family': [],
        'genus': [],
        'species': [],
        'strain': [],
        'sorf': [],
        'slen': [],
        'start-codon': [],
        'protein': [],
        'plen': [],
        'product': [],
        'rbs': [],
        'pfam-hits': [],
        'gravy': [],
        'aromaticity': [],
        'molecular-weight': [],
        'instability': [],
        'isoelectric-point': [],
        'aliphatic-index': [],
        'boman': []
    }

    tar_path: Path = args.output.joinpath(f'sorfdb.{const.SORFDB_VER}.tar.gz')
    csv_path: Path = args.output.joinpath(f'sorfdb.{const.SORFDB_VER}.tsv.gz')
    faa_path: Path = args.output.joinpath(f'sorfdb.{const.SORFDB_VER}.faa.gz')
    fna_path: Path = args.output.joinpath(f'sorfdb.{const.SORFDB_VER}.fna.gz')
    mod_time: float = time()

    with (tarfile.open(tar_path, 'w:gz') as tf,
          xopen(csv_path, 'w', threads=args.threads, compresslevel=9) as tsv_handle,
          xopen(faa_path, 'w', threads=args.threads, compresslevel=9) as faa_handle,
          xopen(fna_path, 'w', threads=args.threads, compresslevel=9) as fna_handle):
        tsv = csv.writer(tsv_handle, delimiter='\t', lineterminator='\n')
        tsv.writerow(list(db_fields.keys()))

        print('Export GenBank sORFs...')
        for datapoint in genbank_extract_db_fields(eio.load_json(args.genbank)):
            tsv.writerow(list(datapoint.values()))
            faa_handle.write(f'>{datapoint["id"]}\n{datapoint["protein"]}\n')
            fna_handle.write(f'>{datapoint["id"]}\n{datapoint["sorf"]}\n')

            datapoint_string = json.dumps(datapoint, indent=0)
            datapoint_io = BytesIO(bytes(datapoint_string, 'utf-8'))
            filename: str = 'sorfdb/' + '.'.join([datapoint['source'], datapoint['accession'], datapoint['protein-id'],
                                                  'json']).replace('..', '.')
            info = tarfile.TarInfo(name=filename)
            info.mtime = mod_time
            info.mode = 0o444
            info.size = len(datapoint_string)
            tf.addfile(tarinfo=info, fileobj=datapoint_io)

        print('Export BLASTX GenBank sORFs...')
        for datapoint in genbank_extract_db_fields(extract_blastx_genbank(args.blastx)):
            tsv.writerow(list(datapoint.values()))
            faa_handle.write(f'>{datapoint["id"]}\n{datapoint["protein"]}\n')
            fna_handle.write(f'>{datapoint["id"]}\n{datapoint["sorf"]}\n')

            datapoint_string = json.dumps(datapoint, indent=0)
            datapoint_io = BytesIO(bytes(datapoint_string, 'utf-8'))
            filename: str = 'sorfdb/' + '.'.join([datapoint['source'], datapoint['accession'], datapoint['protein-id'],
                                                  'json']).replace('..', '.')
            info = tarfile.TarInfo(name=filename)
            info.mtime = mod_time
            info.mode = 0o444
            info.size = len(datapoint_string)
            tf.addfile(tarinfo=info, fileobj=datapoint_io)

        print('Export SmProt sORFs...')
        for datapoint in smprot_extract_db_fields(args.smprot, args.pfam, threads=args.threads):
            # print(datapoint)
            tsv.writerow(list(datapoint.values()))
            faa_handle.write(f'>{datapoint["id"]}\n{datapoint["protein"]}\n')
            if len(datapoint["sorf"]) > 0:
                fna_handle.write(f'>{datapoint["id"]}\n{datapoint["sorf"]}\n')

            datapoint_string = json.dumps(datapoint, indent=0)
            datapoint_io = BytesIO(bytes(datapoint_string, 'utf-8'))
            filename: str = 'sorfdb/' + '.'.join([datapoint['source'], datapoint['accession'], datapoint['protein-id'],
                                                  'json']).replace('..', '.')
            info = tarfile.TarInfo(name=filename)
            info.mtime = mod_time
            info.mode = 0o444
            info.size = len(datapoint_string)
            tf.addfile(tarinfo=info, fileobj=datapoint_io)

        print('Import taxonomy...')
        taxonomy: dict[int, dict[str, str]] = parse_taxon_ids(args.taxonomy)

        for external_db in (args.swissprot, args.uniprot):
            for datapoint in fasta_extract_db_fields(Path(external_db), args.pfam, taxonomy, threads=args.threads):
                tsv.writerow(list(datapoint.values()))
                faa_handle.write(f'>{datapoint["id"]}\n{datapoint["protein"]}\n')

                datapoint_string = json.dumps(datapoint, indent=0)
                datapoint_io = BytesIO(bytes(datapoint_string, 'utf-8'))
                filename: str = 'sorfdb/' + '.'.join([datapoint['source'], datapoint['uid'], datapoint['entry-name'],
                                                      'json']).replace('..', '.')
                info = tarfile.TarInfo(name=filename)
                info.mtime = mod_time
                info.mode = 0o444
                info.size = len(datapoint_string)
                tf.addfile(tarinfo=info, fileobj=datapoint_io)


if __name__ == '__main__':
    main()

# EOF
