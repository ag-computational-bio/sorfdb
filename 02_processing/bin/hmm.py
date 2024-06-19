#!/usr/bin/env python3

from pathlib import Path
from collections import namedtuple
from argparse import ArgumentParser
from typing import Optional, Union

import ijson
import pyhmmer
from xopen import xopen

import extended_io as eio
import constants as const


def get_protein_entries(json_file: Path):
    """
    Yield all protein sequences with their accession number and locus tag.
    :param json_file: Path to sORF json file
    :return: accession, locus, protein
    """
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
                yield accession_number, locus_tag, aa
                locus_tag, aa = '', ''


def get_protein_sequences(json_file: Path) -> list[pyhmmer.easel.DigitalSequence]:
    """
    Extract the sORF aa sequences from the JSON file and convert them to a list of pyhmmer Sequences.
    :param json_file: Path to sORF json file
    :return: list of pyhmmer DigitalSequences from the sORF proteins
    """
    sequences: list[pyhmmer.easel.DigitalSequence] = []
    for accession_number, locus_tag, aa in get_protein_entries(json_file):
        sequences.append(
            pyhmmer.easel.TextSequence(sequence=aa,
                                       accession=bytes(accession_number, 'utf-8'),
                                       name=bytes(locus_tag, 'utf-8')).digitize(pyhmmer.easel.Alphabet.amino())
        )
    return sequences


def perform_hmmsearch(proteins: list[pyhmmer.easel.DigitalSequence],
                      hmm_path: str,
                      cutoff: str,
                      best_hit: bool,
                      threads: int = 1) -> list[namedtuple]:
    """
    Perform a hmmsearch on the given proteins with the HMM profiles in the given file (.h3m)
    :param proteins: list of proteins in pyhmmer DigitalSequence format
    :param hmm_path: path to the .hm3/.hmm file containing the HMM profiles
    :param cutoff: choose cutoff threshold
    :param best_hit: only return best hits
    :param threads: number of threads
    :return: list of best hits (namedtuple)
    """
    # https://pyhmmer.readthedocs.io/en/stable/examples/fetchmgs.html
    Result: namedtuple = namedtuple('Result', ['accession', 'query', 'subject', 'bitscore', 'evalue'])
    results: list[namedtuple] = []

    with pyhmmer.plan7.HMMFile(hmm_path) as hmm:
        for top_hits in pyhmmer.hmmsearch(hmm, proteins, bit_cutoffs=cutoff, cpus=threads):
            for hit in top_hits:
                subject = hit.best_domain.alignment.hmm_name.decode()
                results.append(Result(accession=hit.accession.decode(),
                                      query=hit.name.decode(),
                                      subject=subject,
                                      bitscore=hit.score,
                                      evalue=hit.evalue))
    if best_hit:
        return get_best_hits(results)
    return results


def get_best_hits(results: list[namedtuple]) -> list[namedtuple]:
    """
    Return only the hit with the highest bitscore for each query.
    :param results: list of all hits (namedtuple) of the HMM search
    :return: list of best hits (namedtuple)
    """
    # https://pyhmmer.readthedocs.io/en/stable/examples/fetchmgs.html#Filtering-results
    best_results: dict[str, namedtuple] = {}
    keep_query: set[str] = set()
    for result in results:
        if result.query in best_results:
            previous_bitscore: float = best_results[result.query].bitscore
            if result.bitscore > previous_bitscore:
                best_results[result.query] = result
                keep_query.add(result.query)
            elif result.bitscore == previous_bitscore:
                if best_results[result.query].subject != result.subject:
                    keep_query.remove(result.query)
        else:
            best_results[result.query] = result
            keep_query.add(result.query)
    return [best_results[k] for k in sorted(best_results) if k in keep_query]


def reformat_results(hmmsearch_results: list[namedtuple],
                     prefix: str,
                     evalue: Optional[float] = None) -> dict[str, dict[str, dict[str, list[Union[float, str]]]]]:
    """
    Reformat the hmmsearch results into a dictonary to add to the sORF data.
    :param hmmsearch_results: list of best hits (namedtuple)
    :param prefix: prefix based on the HMM database
    :param evalue: maximum e-value
    :return: {accession_number: {locus_tag: {prefix-evalue: [float, ...], ...}}, ...}
    """
    hmmsearch_hits: dict[str, dict[str, dict[str, list[Union[float, str]]]]] = dict()
    for result in hmmsearch_results:
        if evalue and result.evalue > evalue:
            continue
        hmmsearch_hits.setdefault(result.accession, dict())
        hmmsearch_hits[result.accession].setdefault(result.query, {f"{prefix}-subject": [],
                                                                   f"{prefix}-evalue": []})
        hmmsearch_hits[result.accession][result.query][f"{prefix}-subject"].append(result.subject)
        hmmsearch_hits[result.accession][result.query][f"{prefix}-evalue"].append(result.evalue)
    return hmmsearch_hits


def drop_proteins_with_antifam_hit(proteins: list[pyhmmer.easel.DigitalSequence],
                                   antifam: dict[str, dict[str, dict[str, list[Union[float, str]]]]]) -> list[pyhmmer.easel.DigitalSequence]:
    """
    Drop proteins with AntiFam hits.
    :param proteins: list of pyhmmer DigitalSequences from the sORF proteins
    :param antifam: {accession_number: {locus_tag: {antifam-evalue: [float, ...], ...}}, ...}
    :return: list of pyhmmer DigitalSequences from the sORF proteins without AntiFam hits
    """
    block_indices: set[int] = set()
    for i, protein in enumerate(proteins):
        accession: str = protein.accession.decode('utf-8')
        locus: str = protein.name.decode('utf-8')
        if accession in antifam and locus in antifam[accession]:
            block_indices.add(i)
    return [protein for i, protein in enumerate(proteins) if i not in block_indices]


def drop_antifam_hits(sorf_data: dict,
                      antifam: dict[str, dict[str, dict[str, list[Union[float, str]]]]]) -> dict:
    """
    Drop all sORFs with an Antifam hit.
    :param sorf_data: sORF data
    :param antifam: sORF hits with Antifam
    :return: sORF data without Antifam hits
    """
    i: int = 0
    for accession, antifam_hits in antifam.items():
        if sorf_data[accession]['assembly'] in const.REFERENCE_ACCESSIONS:  # Skip reference genomes
            continue
        for locus in antifam_hits.keys():
            sorf_data[accession]['sorfs'].pop(locus)
            i += 1
    print(f'\tDropped {i} sORFs with AntiFam hits.')
    return sorf_data


def dump_antifam_hits(sorf_data: dict,
                      antifam: dict[str, dict[str, dict[str, list[Union[float, str]]]]]) -> tuple[dict, dict]:
    """
    Drop all sORFs with an Antifam hit.
    :param sorf_data: sORF data
    :param antifam: sORF hits with Antifam
    :return: sORF data without Antifam hits
    """
    antifam_sorfs: dict = dict()
    for accession, antifam_hits in antifam.items():
        if sorf_data[accession]['assembly'] in const.REFERENCE_ACCESSIONS:  # Skip reference genomes
            continue
        for locus in antifam_hits.keys():
            antifam_sorfs.setdefault(accession, {k: v for k, v in sorf_data[accession].items() if k != 'sorfs'})
            antifam_sorfs[accession].setdefault('sorfs', dict())
            antifam_sorfs[accession]['sorfs'][locus] = sorf_data[accession]['sorfs'].pop(locus)
    return sorf_data, antifam_sorfs


def add_hmmsearch_data(data: dict,
                       antifam: dict[str, dict[str, dict[str, list[Union[float, str]]]]],
                       pfam: dict[str, dict[str, dict[str, list[Union[float, str]]]]],
                       skip_antifam: bool = False,
                       inverse: bool = False) -> dict:
    """
    Add the Antifam and Pfam hit data to the sORF data.
    :param data: sORF data from JSON file
    :param antifam: antifam hits data
    :param pfam: pfam hits data
    :param skip_antifam: Skip addibng of Antifam hits (if previously deleted)
    :param inverse: Add Antifam and Pfam hits by traversing the sORF data dictionary. Else: Directly add from the hits.
    :return: enriched sORF data
    """
    if inverse:
        for accession_number, bacteria_data in data.items():
            if 'sorfs' not in bacteria_data:
                continue
            for locus_tag, sorf_data in bacteria_data['sorfs'].items():
                antifam_hit: Optional[dict[str, Union[float, str]]] = get_hit_dict(antifam, accession_number, locus_tag)
                pfam_hit: Optional[dict[str, Union[float, str]]] = get_hit_dict(pfam, accession_number, locus_tag)
                if antifam_hit is not None:
                    sorf_data.update(antifam_hit)
                if pfam_hit is not None:
                    sorf_data.update(pfam_hit)
    else:
        if not skip_antifam:
            for accession, antifam_hits in antifam.items():
                for locus, antifam_hit in antifam_hits.items():
                    try:
                        data[accession]['sorfs'][locus].update(antifam_hit)
                    except KeyError:
                        continue
        for accession, pfam_hits in pfam.items():
            for locus, pfam_hit in pfam_hits.items():
                try:
                    data[accession]['sorfs'][locus].update(pfam_hit)
                except KeyError:
                    continue
    return data


def get_hit_dict(hits: dict[str, dict[str, dict[str, list[Union[float, str]]]]],
                 accession_number: str,
                 locus_tag: str) -> Optional[dict[str, list[Union[float, str]]]]:
    """
    Return the Antifam/Pfam hit data if available
    :param hits: antifam/pfam hits data
    :param accession_number:
    :param locus_tag:
    :return: antifam/pfam hit
    """
    if accession_number in hits and locus_tag in hits[accession_number]:
        return hits[accession_number][locus_tag]


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Extract sORFs from gbff files.')
    parser.add_argument('--data', '-d', type=Path, help='Path to sORF data (JSON file)')
    parser.add_argument('--source', '-s', type=str, help='Data source: annotated, hypothetical or blastx data')
    parser.add_argument('--antifam', '-a', type=Path, help='Path to antifam .h3m file')
    parser.add_argument('--pfam', '-p', type=Path, help='Path to Pfam .h3m file')
    parser.add_argument('--remove_antifam', '-r', action='store_true', help='Remove all proteins with an antifam hit.')
    parser.add_argument('--dump_antifam', '-u', action='store_true', help='Remove and dump all proteins with an antifam hit.')
    parser.add_argument('--output', '-o', type=Path, default=Path('./'), help='Output path')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    return args.data, args.source, args.antifam, args.pfam, args.remove_antifam, args.dump_antifam, args.output, args.threads


def main():
    data_path, source, antifam_path, pfam_path, remove_antifam, dump_antifam, out_path, threads = parse_arguments()
    if remove_antifam and dump_antifam:
        print('Either choose "remove_antifam" OR "dump_antifam".')
        exit(1)

    ijson.get_backend('yajl2_cffi')
    print('Get protein sequences...')
    proteins: list[pyhmmer.easel.DigitalSequence] = get_protein_sequences(data_path)

    print('Start HMM search...')
    print('\tAntiFam')
    antifam: dict[str, dict[str, dict[str, list[Union[float, str]]]]]
    antifam = reformat_results(perform_hmmsearch(proteins,
                                                 antifam_path,
                                                 cutoff='gathering',
                                                 best_hit=True,
                                                 threads=threads),
                               evalue=1E-5,
                               prefix='antifam')

    if remove_antifam:
        proteins = drop_proteins_with_antifam_hit(proteins, antifam)

    print('\tpFam')
    pfam: dict[str, dict[str, dict[str, list[Union[float, str]]]]]
    pfam = reformat_results(perform_hmmsearch(proteins,
                                              pfam_path,
                                              cutoff='gathering',
                                              best_hit=False,
                                              threads=threads),
                            evalue=1E-5,
                            prefix='pfam')

    del proteins

    if remove_antifam:
        print('Load sORFs data & drop Antifam sORFs...')
        data = drop_antifam_hits(sorf_data=eio.load_json(data_path),
                                 antifam=antifam)
    elif dump_antifam:
        print('Load sORFs data & dump Antifam sORFs...')
        data, antifam_sorf = dump_antifam_hits(sorf_data=eio.load_json(data_path),
                                               antifam=antifam)
        antifam_sorf_out_path: Path = out_path.joinpath(f"antifam.{source}.json.gz")
        eio.export_to_json(add_hmmsearch_data(antifam_sorf,
                                              antifam=antifam,
                                              pfam=pfam),
                           antifam_sorf_out_path,
                           threads=threads)
        del antifam_sorf
    else:
        print('Load sORFs data...')
        data: dict = eio.load_json(data_path)

    print('Write out...')
    data_out_path: Path = out_path.joinpath(f"{source}.json")
    eio.export_to_json(add_hmmsearch_data(data,
                                          antifam=antifam,
                                          pfam=pfam,
                                          skip_antifam=(remove_antifam or dump_antifam)),
                       data_out_path,
                       threads=threads)


if __name__ == '__main__':
    main()

# EOF
