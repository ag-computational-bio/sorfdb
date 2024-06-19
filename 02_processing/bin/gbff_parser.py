#!/usr/bin/env python3

import os
import re
import subprocess as sp
from pathlib import Path
from copy import deepcopy
from argparse import ArgumentParser
from typing import Any, Optional, Union

import pyrodigal
from xopen import xopen
from Bio import Seq, SeqIO

import extended_io as eio
import biofunctions as bf
import constants as const

NUCLEOTIDES: set[str] = {'A', 'C', 'G', 'T'}


def gbff_parser(gbff_path: Path,
                fasta_path: Path,
                assembly: str,
                threads: int = 1) -> tuple[dict[str, dict],
                                           dict[str, dict],
                                           dict[str, Union[int, dict[str, Any]]],
                                           pyrodigal.OrfFinder]:
    """
    Extract all complete and hypothetical sORF entries
    from genbank ('genomes/genbank/bacteria/') *.gbff.gz files.
    :param gbff_path: Path to *.gbff.gz file
    :param fasta_path: Output path for the fasta file
    :param assembly: Genome assembly number
    :param threads: Number of threads to use
    :return: Dictionaries with extracted valid/hypothetical sORFs for each record and feature ranges.
    """
    valid_sorfs: dict[str, dict] = dict()
    hypothetical_sorfs: dict[str, dict] = dict()

    annotations: dict[str, Union[int, dict[str, Any]]] = dict()
    translation_tables: dict[int, int] = dict()

    fna_output_handle = xopen(fasta_path, mode='w', compresslevel=9, threads=threads)

    orf_finder: pyrodigal.OrfFinder = initialize_pyrodigal(gbff_path)
    reference_genome: bool = is_reference_genome(gbff_path)
    representative_genome: bool = is_representative_genome(gbff_path)
    pgap_annotated: bool = is_pgap_annotated(gbff_path)

    with xopen(gbff_path, mode='r') as f:
        for record in SeqIO.parse(f, 'genbank'):
            SeqIO.write(record, fna_output_handle, 'fasta')

            topology: str = record.annotations['topology']
            taxonomy: tuple[str, tuple[str, ...]] = get_taxonomy(record.annotations)
            protein_id_duplicate_counter: int = 1

            annotations.setdefault(record.id, {'features': [],
                                               'genome': str(record.seq).upper(),
                                               'topology': topology,
                                               'taxonomy': taxonomy})

            base_entry: dict = {'taxonomy': taxonomy,
                                'assembly': assembly,
                                'sorfs': dict()}

            for feat in record.features:
                strand: int = feat.location.strand

                if feat.type != 'source':
                    annotations[record.id]['features'].append({'start': int(feat.location.start),
                                                               'stop': int(feat.location.end),
                                                               'strand': strand,
                                                               'type': feat.type})

                # Flags: https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
                # sORF
                if feat.type == 'CDS' and int(feat.location.end) - int(feat.location.start) <= const.MAX_SORF_LENGTH and \
                        is_triplet_sequence(feat) and not is_fuzzy(feat) and not is_pseudogene(feat) and not is_fragment(feat):

                    product: str = feat.qualifiers['product'][0] if 'product' in feat.qualifiers else None
                    if 'protein_id' in feat.qualifiers:
                        tag: str = feat.qualifiers['protein_id'][0]
                    elif 'locus_tag' in feat.qualifiers:
                        tag: str = feat.qualifiers['locus_tag'][0]
                    else:
                        tag: str = f"{feat.type}_{feat.location.start}_{feat.location.end}_{feat.location.strand}"

                    sorf: str = feat.extract(str(record.seq.upper()))
                    if not bf.is_dna(sorf) or not has_stop_codon(sorf):
                        continue

                    upstream, downstream = extract_flanking_sequences(str(record.seq.upper()),
                                                                      strand,
                                                                      int(feat.location.start),
                                                                      int(feat.location.end),
                                                                      topology)

                    if 'transl_table' in feat.qualifiers:
                        transl_table: int = int(feat.qualifiers['transl_table'][0])
                    elif 'transl_table' not in locals():
                        transl_table: Optional[int] = identify_translation_table(sorf)
                    if transl_table not in const.TRANSLATION_TABLES:
                        del transl_table
                        continue

                    translation_tables.setdefault(transl_table, 0)
                    translation_tables[transl_table] += 1

                    rbs: Optional[dict] = predict_rbs(sorf, upstream, downstream, orf_finder)

                    if 'translation' in feat.qualifiers:
                        aa: str = feat.qualifiers['translation'][0]
                    else:
                        aa: str = f"M{Seq.translate(sorf, table=transl_table, stop_symbol='')[1:]}"

                    tag = replace_invalid_characters(tag)

                    if product is not None and len(product.strip()) > 1:
                        pl = product.lower()
                        hypothetical_product: bool = is_hypothetical_protein(pl)
                    else:
                        product = 'hypothetical protein'
                        hypothetical_product: bool = True

                    sorf_entry: dict[str, Union[str, int]] = {'transl_table': transl_table,
                                                              'product': product,
                                                              'orf': sorf,
                                                              'up': upstream,
                                                              'down': downstream,
                                                              'aa': aa}
                    if rbs is not None:
                        sorf_entry.update(rbs)

                    if (reference_genome or representative_genome or pgap_annotated) and not hypothetical_product:
                        if record.id in valid_sorfs and tag in valid_sorfs[record.id]['sorfs']:
                            tag = f'{tag}#{protein_id_duplicate_counter}'
                            protein_id_duplicate_counter += 1

                        valid_sorfs.setdefault(record.id, deepcopy(base_entry))
                        valid_sorfs[record.id]['sorfs'][tag] = sorf_entry
                        valid_sorfs[record.id]['sorfs'][tag]['type'] = 'annotated'
                    else:
                        if record.id in hypothetical_sorfs and tag in hypothetical_sorfs[record.id]['sorfs']:
                            tag = f'{tag}#{protein_id_duplicate_counter}'
                            protein_id_duplicate_counter += 1

                        hypothetical_sorfs.setdefault(record.id, deepcopy(base_entry))
                        hypothetical_sorfs[record.id]['sorfs'][tag] = sorf_entry
                        hypothetical_sorfs[record.id]['sorfs'][tag]['type'] = 'hypothetical'

    fna_output_handle.close()

    if len(annotations) > 0:
        annotations['transl_table'] = most_used_transl_table(translation_tables)

    print(f'\tAnnotated sORFs: {sum([len(bacteria_data["sorfs"]) for bacteria_data in valid_sorfs.values()])}')
    print(f'\tHypothetical sORFs: {sum([len(bacteria_data["sorfs"]) for bacteria_data in hypothetical_sorfs.values()])}')
    return valid_sorfs, hypothetical_sorfs, annotations, orf_finder


def extract_flanking_sequences(sequence: str,
                               strand: int,
                               start: int,
                               stop: int,
                               topology: str,
                               sorf: str = None,
                               upstream_length: int = const.UPSTREAM_LENGTH,
                               downstream_length: int = const.DOWNSTREAM_LENGTH,
                               trim_degenerated: bool = False) -> tuple[str, str]:
    """
    Extract the sORF ORF and up/downstream sequences
    :param sequence: DNA sequence of the genome
    :param strand: DNA strand of the ORF
    :param start: Start position of the ORF
    :param stop: Stop position of the ORF
    :param topology: Genome topology
    :param sorf: sORF sequence for assertion test
    :param upstream_length: target length for the upstream sequence
    :param downstream_length: target length for the downstream sequence
    :param trim_degenerated: Leading and trailing degenerated nucleotides
    :return: Upstream, ORF, downstream sequences
    """
    # Extract upstream sequence
    up_length: int = upstream_length if strand == 1 else downstream_length

    if start >= up_length:
        upstream: str = sequence[start - up_length:start]
        up_position = start - up_length
    elif topology == 'circular':
        upstream: str = sequence[len(sequence) - (up_length - start):] + sequence[:start]
        up_position = len(sequence) - (up_length - start)
    else:
        upstream: str = sequence[:start]
        up_position = 0

    # Extract downstream sequence
    down_length: int = downstream_length if strand == 1 else upstream_length

    if stop <= (len(sequence) - down_length):
        downstream: str = sequence[stop:stop + down_length]
        down_position = stop + down_length
    elif topology == 'circular':
        downstream: str = sequence[stop:] + sequence[:down_length - (len(sequence) - stop)]
        down_position = down_length - (len(sequence) - stop)
    else:
        downstream: str = sequence[stop:]
        down_position = len(sequence)

    # Create complement if needed
    if strand == -1:
        upstream, downstream = Seq.reverse_complement(downstream), Seq.reverse_complement(upstream)
        if sorf is not None and up_position < start and down_position > stop:
            assert upstream + sorf + downstream == Seq.reverse_complement(sequence[up_position:down_position])
    if trim_degenerated:
        # Strip leading degenerate bases in the upstream sequence and trailing degenerate bases in the downstream sequence
        degenerated: str = 'NWSMKRYBDHV'
        upstream, downstream = upstream.lstrip(degenerated), downstream.rstrip(degenerated)
        return ltrim(upstream), rtrim(downstream)
    return upstream, downstream


def ltrim(seq: str) -> str:
    """
    Trim all nucleotides on the left in a sequence, starting from the ambiguous nucleotide furthest to the right.
    :param seq: DNA sequence; can contain ambiguous nucleotides
    :return: DNA sequence without ambiguous nucleotides
    """
    for i in range(len(seq)-1, 0, -1):
        if seq[i] not in NUCLEOTIDES:
            return seq[i+1:]
    return seq


def rtrim(seq: str) -> str:
    """
    Trim all nucleotides on the right in a sequence, starting from the ambiguous nucleotide furthest to the left.
    :param seq: DNA sequence; can contain ambiguous nucleotides
    :return: DNA sequence without ambiguous nucleotides
    """
    for i, nucleotide in enumerate(seq):
        if nucleotide not in NUCLEOTIDES:
            return seq[:i]
    return seq


def initialize_pyrodigal(gbff_path: Path) -> pyrodigal.OrfFinder:
    """
    Create a pyrodigal.OrfFinder instance trained on all contigs if contigs >= 100000bp else use meta mode.
    :param gbff_path: path to genbank file
    :return: pyrodigal.OrfFinder instance
    """
    min_contig_length: int = 100000
    transl_table: int = 11
    sequences = []
    contig_lenth: int = 0

    with xopen(gbff_path, mode='r') as f:
        for i, record in enumerate(SeqIO.parse(f, 'genbank')):
            if i == 0:
                transl_table = peek_translation_table(record.features)
            sequences.append(str(record.seq))
            contig_lenth += len(record.seq)

    if contig_lenth >= min_contig_length:
        orf_finder = pyrodigal.OrfFinder(min_gene=const.MIN_SORF_LENGTH,
                                         max_overlap=const.MIN_SORF_LENGTH,
                                         closed=True)
        orf_finder.train(*sequences, translation_table=transl_table)
    else:
        orf_finder = pyrodigal.OrfFinder(min_gene=const.MIN_SORF_LENGTH,
                                         max_overlap=const.MIN_SORF_LENGTH,
                                         meta=True,
                                         closed=True)
    return orf_finder


def shorten_upstream_sequence(seq: str, length: int) -> str:
    """
    Shorten the upstream nucleotide sequence to the given length.
    :param seq: Upstream DNA sequence
    :param length: Maximum sequence length
    :return: Shortened upstream sequence
    """
    if len(seq) > length:
        return seq[-length:]
    return seq


def shorten_downstream_sequence(seq: str, length: int) -> str:
    """
    Shorten the downstream nucleotide sequence to the given length.
    :param seq: Downstream DNA sequence
    :param length: Maximum sequence length
    :return: Shortened downstream sequence
    """
    if len(seq) > length:
        return seq[:max(0, length)]
    return seq


def predict_rbs(sorf: str, upstream: str, downstream: str,
                orf_finder: pyrodigal.OrfFinder) -> Optional[dict[str, Union[str, int]]]:
    """
    Predict RBS of sORFs using Pyrodigal.
    :param sorf: sORF sequence
    :param upstream: upstream sequence
    :param downstream: downstream sequence
    :param orf_finder: orf_finder instance of pyrodigal
    :return: RBS data if existing
    """
    up: str = shorten_upstream_sequence(upstream, const.RBS_UPSTREAM_LENGTH)
    down: str = shorten_downstream_sequence(downstream, 50)
    sorf_start: int = len(up) + 1
    sorf_stop: int = len(up) + len(sorf)
    rerun: int = 3  # allow 2 reruns; for reversed strand hits

    while rerun > 0:
        # shortest RBS = 6b (3bp SD + 3bp spacer)
        if len(up) < 6:
            return None

        # try to extract the gene
        try:
            genes: pyrodigal.Genes = orf_finder.find_genes(f"{up}{sorf}{down}")
        except RuntimeError:
            print('Pyrodigal: RBS prediction failed')
            return None
        if len(genes) == 0:
            break
        for prediction in genes:
            if prediction.strand == 1 and prediction.start_type != 'Edge':
                if prediction.begin == sorf_start and prediction.end == sorf_stop:
                    return {'uses_sd': 1 if prediction.rbs_motif is not None else 0,
                            'rbs_motif': prediction.rbs_motif,
                            'rbs_spacer': prediction.rbs_spacer}
            if len(genes) == 1:
                if prediction.strand == 1:
                    if prediction.begin < sorf_start:
                        up = shorten_upstream_sequence(up, len(up) - prediction.begin)
                        sorf_start = len(up) + 1
                        sorf_stop = len(up) + len(sorf)
                    if prediction.end > sorf_stop:
                        down = shorten_downstream_sequence(down, prediction.end - sorf_stop)
                else:
                    if rerun - 2 > 0:
                        cutoff: int = len(down) - len(f"{up}{sorf}{down}") - prediction.end + 1
                        down = shorten_downstream_sequence(down, cutoff)
                    else:
                        down = ''
        rerun -= 1


def is_fragment(feat) -> bool:
    """
    Check if a CDS encodes a protein fragment or truncated protein.
    :param feat: Bio.SeqFeature
    :return: Boolean
    """
    if 'product' in feat.qualifiers:
        product = feat.qualifiers['product'][0].lower()
        return 'fragment' in product or 'truncat' in product or 'partial' in product
    return False


def is_fuzzy(feat) -> bool:
    """
    Check if a feature has a fuzzy position or is a join or ordered sequence.
    :param feat: Bio.SeqFeature; contains location
    :return: Boolean
    """
    location = str(feat.location)
    return '<' in location or '>' in location or location.startswith('join') or location.startswith('order')


def is_pseudogene(feat) -> bool:
    """
    Check if a feature is a pseudogene.
    :param feat: Bio.SeqFeature
    :return: Boolean
    """
    if feat.type == 'CDS':
        return 'pseudo' in feat.qualifiers or 'pseudogene' in feat.qualifiers or \
               ('product' in feat.qualifiers and 'pseudogene' in feat.qualifiers['product'][0])
    return False


def is_triplet_sequence(feat) -> bool:
    """
    Check if a sequence is a valid codon triplet sequence (is a multiple of 3)
    :param feat: Bio.SeqFeature; contains location
    :return: Boolean
    """
    return (int(feat.location.end) - int(feat.location.start)) % 3 == 0


def has_stop_codon(seq: str) -> bool:
    """
    Check if the sequence ends with a stop codon (TAG, TAA, TGA)
    :param seq: DNA sequence
    :return: Has a stop codon
    """
    return seq.endswith('TAG') or seq.endswith('TAA') or seq.endswith('TGA')


def is_hypothetical_protein(product: str) -> bool:
    """
    Check if a protein product is hypothetical.
    :param product: product annotation
    :return: is hypothetical
    """
    if 'hypothetical' in product or 'putative' in product or 'unknown' in product or 'possible' in product or \
            'uncharacterized' in product or 'uncharacterised' in product or 'probable' in product or \
            'dubious' in product or 'doubtful' in product or 'questionable' in product or product == 'protein':
        return True
    return False


def identify_translation_table(orf: str) -> Optional[int]:
    """
    Try to determine the correct translation table for a bacterial DNA sequence.
    Cannot discern table 4 and 25.
    :param orf: Bacterial DNA sequence
    :return: Optional: Table number
    """
    if orf[-3:] in const.STOP_CODONS:  # check ending stop codon
        for table in (11, 4):
            aa: str = f"M{Seq.translate(orf, table=table, stop_symbol='')[1:]}"
            if len(aa) == int((len(orf)-3)/3):  # check internal stop codons
                return table


def replace_invalid_characters(s: str) -> str:
    """
    Replace whitespace characters with an underscore '_'.
    :param s: input string
    :return: string with whitespaces replaced by '_'
    """
    return re.sub(r"(\s+|%[0-9][A-Z]\s?)", '_', s)


def most_used_transl_table(translation_tables: dict[int, int]) -> int:
    """
    Return the most used translation table from the dictionary {table: counter}.
    :param translation_tables:  Dictionary {table: counter}
    :return: Translation table
    """
    translation_table: int = 11  # default
    max_count: int = 0
    for table, count in translation_tables.items():
        try:
            table = int(table)
        except TypeError:
            continue

        if count > max_count:
            max_count = count
            translation_table = table
    if translation_table not in const.TRANSLATION_TABLES:
        raise ValueError(f'Invalid translation table: {translation_table}')
    return translation_table


def is_pgap_annotated(gbff_path: Path) -> bool:
    """
    Check if the input genome is PGAP annotated.
    :param gbff_path: Path to the genome
    :return: is PGAP annotated
    """
    with xopen(gbff_path, mode='r') as f:
        for record in SeqIO.parse(f, 'genbank'):
            if 'comment' in record.annotations and \
                    ('prokaryotic genome annotation pipeline' in record.annotations['comment'].lower().replace('\n', ' ') or
                     'pgap' in record.annotations['comment'].lower().replace('\n', ' ')):
                return True
            if 'structured_comment' in record.annotations and 'Genome-Annotation-Data' in record.annotations['structured_comment']:
                if 'Annotation Pipeline' in record.annotations['structured_comment']['Genome-Annotation-Data'] and \
                        'prokaryotic genome annotation pipeline' in record.annotations['structured_comment']['Genome-Annotation-Data']['Annotation Pipeline'].lower().replace('\n', ' '):
                    return True
                if 'Annotation Method' in record.annotations['structured_comment']['Genome-Annotation-Data'] and \
                        ('reference protein' in record.annotations['structured_comment']['Genome-Annotation-Data']['Annotation Method'].lower().replace('\n', ' ')
                         or 'genemarks+' in record.annotations['structured_comment']['Genome-Annotation-Data']['Annotation Method'].lower().replace('\n', ' ')):
                    return True
    return False


def is_reference_genome(gbff_path: Path) -> bool:
    """
    Check if the input genome is a reference genome.
    :param gbff_path: Path to the genome
    :return: Is reference genome
    """
    return str(gbff_path).split('/')[-1].split('.')[0] in const.REFERENCE_ACCESSIONS


def is_representative_genome(gbff_path: Path) -> bool:
    """
    Check if the input genome is a representative genome.
    :param gbff_path: Path to the genome
    :return: Is representative genome
    """
    return str(gbff_path).split('/')[-1].split('.')[0] in const.REPRESENTATIVE_ACCESSIONS


def get_taxonomy(annotations: dict) -> tuple[str, tuple[str, ...]]:
    """
    Get the taxonomic entry of a genbank record.
    :param annotations: genbank record annotations
    :return: taxonomy: ('Species', ('Bacteria', ...))
    """
    organism: str = annotations['organism']
    tree: list[str] = [level.capitalize() for level in annotations['taxonomy']]
    if len(tree) == 0:
        tree = ['Bacteria', 'Unknown']
    elif tree[0] != 'Bacteria':
        if 'Bacteria' in tree:
            i: int = tree.index('Bacteria')
            tree = tree[i:]
        else:
            tree.insert(0, 'Bacteria')
    return organism, tuple(tree)


def peek_translation_table(record_features: list) -> int:
    """
    Peek the first translation table of a genbank record.
    :param record_features: Features of the accession sequence
    :return: first translation table found
    """
    for feat in record_features:
        if feat.type == 'CDS':
            if 'transl_table' in feat.qualifiers:
                return int(feat.qualifiers['transl_table'][0])
    return 11


def is_annotated(gbff_path: Path) -> bool:
    """
    Check if a genbank genome provides any CDS annotations.
    :param gbff_path: path to genbank file
    :return: Contains annotated CDS
    """
    with xopen(gbff_path, mode='r') as f:
        for record in SeqIO.parse(f, 'genbank'):
            if len(record.id) == 0:
                ValueError('Accession number (record.id) must be given for a genbank entry.')
            for feat in record.features:
                if feat.type == 'CDS':
                    return True
    return False


def blastx(dmnd_db_path: Path, fna_path: Path, assembly: str, out_path: Path,
           annotations: dict[str, Union[int, dict[str, Any]]],
           orf_finder: pyrodigal.OrfFinder, threads: int = 1) -> dict[str, dict]:
    """
    Use Daimond blastx to find missed sORFs in the annotated genomes. Check them to encode a valid protein sequence and
    add them to the hypothetical dataset to verify their Score Ratio values.
    :param dmnd_db_path: Path to diamond database with known sORFs
    :param fna_path: Path to FNA file of the used genome (extracted from the genbank file)
    :param assembly: Genome assembly number
    :param out_path: Output path
    :param annotations: Dictionary containing all annotations and additional data from the genbank file
    :param orf_finder: pyrodigal.OrfFinder instance
    :param threads: Number of threads
    :return: Dictionary containing extracted hypothetical sORFs and missed hypothetical sORFs
    """
    blastx_sorfs: dict[str, dict] = dict()
    diamond_out_path: Path = out_path.joinpath(f"{assembly}.blastx.tsv.gz")

    cmd = [
        'diamond',
        'blastx',
        '--db', str(dmnd_db_path),
        '--query', str(fna_path),
        '--out', str(diamond_out_path),
        '--query-gencode', str(annotations['transl_table']),
        '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 'bitscore', 'full_sseq', 'qseq_translated',
        '--matrix', 'BLOSUM62',
        '--ultra-sensitive',
        '--long-reads',
        '--min-orf', '1',
        '--masking', '0',  # -seg
        '--comp-based-stats', '0',  # -composition_base
        '--evalue', '10',
        '--threads', str(threads),
        '--index-chunks', '1',
        '--compress', '1'
    ]

    # print('Run Diamond BLASTX:\n' + ' '.join(cmd))
    proc = sp.run(
        cmd,
        env=os.environ.copy(),
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if proc.returncode != 0:
        print('Command:\n' + ' '.join(cmd))
        print('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        raise Exception(f'Diamond error! error code: {proc.returncode}')

    valid_count: int = 0
    hit_count = 0
    unique_hits: set[str] = set()
    with xopen(diamond_out_path, mode='r') as fh:
        for line in fh:
            (accession_number, cluster_id, identity, alignment_length, query_start, query_end,
             bitscore, cluster_sequence, qseq_translated) = line.strip().split('\t')

            query_start, query_end = int(query_start), int(query_end)
            if query_start < query_end:
                strand: int = 1
                query_start, query_end = query_start - 1, query_end + 3  # -1 for 0 based index, +3 for stop codon
            else:
                strand: int = -1
                query_start, query_end = query_end - 1 - 3, query_start  # -1 for 0 based index, -3 for stop codon

            if not overlaps_features(accession_number, query_start, query_end, strand, annotations):
                hit_count += 1

                sorf: str = annotations[accession_number]['genome'][query_start:query_end]
                if strand == -1:
                    sorf = Seq.reverse_complement(sorf)

                qseq_translated = f"M{qseq_translated[1:]}"

                hit_id: str = f"{accession_number}|{query_start}|{query_end}|{strand}"
                # if not unique or not triplet or has gaps or not sprot length or not acgt or has internal stop codon

                if hit_id in unique_hits or not len(sorf) % 3 == 0 or (query_end - query_start - 3) // 3 != len(qseq_translated) \
                        or not const.MIN_SORF_LENGTH <= len(sorf) <= const.MAX_SORF_LENGTH or sorf[:3] not in const.START_CODONS or \
                        not bf.is_dna(sorf) or '*' in qseq_translated:
                    unique_hits.add(hit_id)
                    continue
                unique_hits.add(hit_id)

                aa_with_stop: str = f"M{Seq.translate(sorf, table=annotations['transl_table'], stop_symbol='*')[1:]}"
                aa: str = f"M{Seq.translate(sorf, table=annotations['transl_table'], stop_symbol='')[1:]}"

                if aa_with_stop.endswith('*') and aa_with_stop.count('*') == 1 and aa == qseq_translated:
                    upstream, downstream = extract_flanking_sequences(annotations[accession_number]['genome'],
                                                                      strand,
                                                                      query_start,
                                                                      query_end,
                                                                      annotations[accession_number]['topology'])

                    rbs: Optional[dict] = predict_rbs(sorf, upstream, downstream, orf_finder)

                    sorf_entry: dict[str, Union[str, int]] = {'strand': strand,
                                                              'start': query_start,
                                                              'stop': query_end,
                                                              'transl_table': annotations['transl_table'],
                                                              'product': 'hypothetical protein',
                                                              'orf': sorf,
                                                              'up': upstream,
                                                              'down': downstream,
                                                              'aa': aa,
                                                              'type': 'blastx'}

                    if rbs is not None:
                        sorf_entry.update(rbs)

                    tag: str = f"sORF_{query_start}_{query_end}_{'+' if strand == 1 else '-'}"
                    blastx_sorfs.setdefault(accession_number, {'taxonomy': annotations[accession_number]['taxonomy'],
                                                               'assembly': assembly,
                                                               'sorfs': dict()})
                    blastx_sorfs[accession_number]['sorfs'][tag] = sorf_entry

                    valid_count += 1
    print('\tNon-overlapping hits:', hit_count, '\n\tValid hits:', valid_count)
    return blastx_sorfs


def overlaps_features(accession_number: str, start: int, stop: int, strand: int,
                      features: dict[str, Union[int, dict[str, Any]]]) -> bool:
    """
    Check if a hit overlaps with an existing feature.
    :param accession_number: Accession number (read) of the hit
    :param start: Start position
    :param stop: Stop position
    :param strand: Strand direction
    :param features: Dictionary containing all features of a genbank entry
    :return: Boolean
    """
    for feat in features[accession_number]['features']:
        hit_frame: int = start % 3
        feat_frame: int = feat['start'] % 3
        if feat['type'] == 'CDS':
            if strand == feat['strand']:
                if hit_frame == feat_frame:
                    # check partial overlap
                    if not (feat['stop'] < start or feat['start'] > stop):
                        return True
                elif start >= feat['start'] and stop <= feat['stop']:
                    # out-of-frame sorf completely overlapped by CDS
                    return True
            else:
                if hit_frame == feat_frame and not (feat['stop'] < start or feat['start'] > stop):
                    return True
                elif start >= feat['start'] and stop <= feat['stop']:
                    # out-of-frame sorf completely overlapped by CDS
                    return True
        else:
            if not (feat['stop'] < start or feat['start'] > stop):
                return True
    return False


def get_assembly_name(gbff_file: Path) -> str:
    """
    Extract the assembly name/number from the file name.
    :param gbff_file: path to annotated genbank file
    :return: assembly name
    """
    assembly_tmp: list[str] = str(gbff_file).split('/')[-1].split('.')
    if assembly_tmp[1].endswith('genomic') and assembly_tmp[1].split('_')[0].isdigit():
        return '.'.join([assembly_tmp[0], assembly_tmp[1].split('_')[0]])
    return assembly_tmp[0]


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Extract sORFs from gbff files.')
    parser.add_argument('--id', '-i', type=str, help='ID of the given file')
    parser.add_argument('--gbff', '-g', type=Path, help='GBFF Genbank file')
    parser.add_argument('--db', '-d', type=Path, help='sORF database (.dmnd)')
    parser.add_argument('--output', '-o', type=Path, default=Path('./'), help='Output path')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    return args.id, args.gbff, args.db, args.output, args.threads


def main():
    file_id, gbff_file, db_path, out_path, threads = parse_arguments()
    fna_path: Path = out_path.joinpath(f"{file_id}.fna.gz")

    print(f"ID: {file_id}\nGBFF: {gbff_file}\nDB: {db_path}\nOUT: {out_path}\nTHREADS: {threads}")
    if not is_annotated(gbff_path=gbff_file):
        print('Not annotated. Exiting.')
        exit(0)

    assembly: str = get_assembly_name(gbff_file)

    print('Extracting sORFs...')
    valid_sorfs, hypothetical_sorfs, annotations, orf_finder = gbff_parser(gbff_file,
                                                                           fasta_path=fna_path,
                                                                           assembly=assembly,
                                                                           threads=threads)

    print('Searching for additional sORFs...')
    blastx_sorfs: dict[str, dict] = dict()
    if len(annotations) > 0:
        blastx_sorfs = blastx(dmnd_db_path=db_path,
                              fna_path=fna_path,
                              assembly=assembly,
                              out_path=out_path,
                              annotations=annotations,
                              orf_finder=orf_finder,
                              threads=threads)

    print('Write out found sORFs...')
    if len(valid_sorfs) > 0:
        eio.export_to_json(valid_sorfs,
                           out_path.joinpath(f"{file_id}.annotated.json.gz"),
                           single_line=False,
                           threads=threads)

    if len(hypothetical_sorfs) > 0:
        eio.export_to_json(hypothetical_sorfs,
                           out_path.joinpath(f"{file_id}.hypothetical.json.gz"),
                           single_line=False,
                           threads=threads)

    if len(blastx_sorfs) > 0:
        eio.export_to_json(blastx_sorfs,
                           out_path.joinpath(f"{file_id}.blastx.json.gz"),
                           single_line=False,
                           threads=threads)
    print('Done')


if __name__ == '__main__':
    main()

# EOF
