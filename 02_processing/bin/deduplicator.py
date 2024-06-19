#!/usr/bin/env python3

import sys
from pathlib import Path
from datetime import datetime

from xopen import xopen

import extended_io as io


data: dict[str, dict] = dict()


def collapse_multihits(mode: str) -> int:
    global data
    duplicates: dict[str, list[tuple[str, str]]] = dict()
    for accession_number, bacteria_data in data.items():
        if 'sorfs' not in bacteria_data:
            continue
        for locus_tag, sorf_data in bacteria_data['sorfs'].items():
            if sorf_data['orf'] not in duplicates:
                duplicates[sorf_data['orf']] = [(accession_number, locus_tag)]
            else:
                first_accession_number, first_tag = duplicates[sorf_data['orf']][0]
                first_sorf_data: dict = data[first_accession_number]['sorfs'][first_tag]
                len_first_complete_seq: int = len(first_sorf_data['up']) + len(first_sorf_data['orf']) + len(first_sorf_data['down'])
                len_new_complete_seq: int = len(sorf_data['up']) + len(sorf_data['orf']) + len(sorf_data['down'])

                if len_first_complete_seq < len_new_complete_seq:
                    duplicates[sorf_data['orf']].insert(0, (accession_number, locus_tag))
                    data[accession_number]['sorfs'][locus_tag]['duplicate'] = duplicates[sorf_data['orf']][1:]
                    # mark as type=annotated if a duplicate is type='annotated'
                    if data[accession_number]['sorfs'][locus_tag]['type'] == 'hypothetical':
                        for duplicate in duplicates[sorf_data['orf']][1:]:
                            duplicate_accession_number, duplicate_tag = duplicate
                            if data[duplicate_accession_number]['sorfs'][duplicate_tag]['type'] == 'annotated':
                                data[accession_number]['sorfs'][locus_tag]['type'] = 'annotated'
                                break
                else:
                    duplicates[sorf_data['orf']].append((accession_number, locus_tag))
                    data[first_accession_number]['sorfs'][first_tag].setdefault('duplicate', [])
                    data[first_accession_number]['sorfs'][first_tag]['duplicate'].append((accession_number, locus_tag))
                    # mark as type=annotated if a duplicate is type='annotated'
                    if data[first_accession_number]['sorfs'][first_tag]['type'] == 'hypothetical':
                        if data[accession_number]['sorfs'][locus_tag]['type'] == 'annotated':
                            data[first_accession_number]['sorfs'][first_tag]['type'] = 'annotated'

    i: int = 0
    for deletables in duplicates.values():
        if len(deletables) > 1:
            i += len(deletables) - 1  # -1 for first valid entry
            for accession_number, tag in deletables[1:]:
                data[accession_number]['sorfs'].pop(tag)
    return i


def total_number_of_sorfs(sorfs: dict[str, dict]) -> int:
    """
    Return the total number of sORFs entries in the given database dictionary.
    :param sorfs: sORF dictionary
    :return: Number of sORF
    """
    i: int = 0
    for sorf_ls in sorfs.values():
        if 'sorfs' in sorf_ls:
            i += len(sorf_ls['sorfs'])
    return i


def export_faa(file_path: str, threads: int):
    """
    Extract all aa sequences from the DB JSON file(s) write them to a compressed FASTA file.
    :param file_path: Path to the output file
    :param threads: Number of threads to use for compression
    """
    global data

    with xopen(file_path, mode='w', compresslevel=9, threads=threads) as f:
        for accession_number, entry in data.items():
            if accession_number == 'metadata':
                continue
            for locus_tag, sorf_data in entry['sorfs'].items():
                header: str = f"{accession_number}|{locus_tag}"
                f.write(f">{header}\n{sorf_data['aa']}\n")


def main():
    _, in_path, mode, source, threads = sys.argv
    in_path = Path(in_path)
    out_path = Path('./')
    threads: int = int(threads)

    if mode not in {'valid', 'hypothetical'}:
        raise ValueError(f"mode={mode} is not valid. Use 'valid' or 'hypothetical'.")

    print('Import data...')
    global data
    for i, json_file in enumerate(io.find_files(in_path, r'.*\.json.*')):
        print(json_file)
        if i == 0:
            data = io.load_json(json_file)
        else:
            for accession_number, bacteria_data in io.load_json(json_file).items():
                if accession_number not in data:
                    data[accession_number] = bacteria_data
                else:
                    for locus_tag, sorf_data in bacteria_data['sorfs'].items():
                        data[accession_number]['sorfs'][locus_tag] = sorf_data

    total: int = total_number_of_sorfs(data)

    if mode == 'valid':
        print('Export all data...')
        all_export_path: str = str(out_path.joinpath("all.valid.json"))
        io.export_to_json(data, all_export_path, threads=threads, compresslevel=9)

    print('Deduplicate data...')
    num_deleted: int = collapse_multihits(mode)

    total_deduplicated: int = total_number_of_sorfs(data)
    data['metadata'] = {'data': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        'total deduplicated': total_deduplicated,
                        'duplicates (dropped)': num_deleted}

    print('Found {} sORFs total.'.format(total))
    print('Found {} unique sORFs.'.format(data['metadata']['total deduplicated']))
    print('Found {} double entries.'.format(num_deleted))

    print('Export deduplicated data...')
    export_path: str = str(out_path.joinpath(f"deduplicated.{mode}.{source}.json"))
    io.export_to_json(data, export_path, threads=threads)


if __name__ == '__main__':
    main()

# EOF
