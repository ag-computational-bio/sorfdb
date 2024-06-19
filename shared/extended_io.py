
import os
import json
from pathlib import Path
from typing import Generator, Optional, Union

from xopen import xopen


def load_json(file_path: Union[str, Path]) -> dict:
    """
    Load a (compressed) JSON file into a dictionary.
    :param file_path: Path to (compressed) JSON file
    :return: Dictionary from given file
    """
    with xopen(file_path, mode='r') as f:
        return json.load(f)


def export_to_json(data: dict, file_path: Union[str, Path], single_line: bool = False, threads: Optional[int] = None,
                   compresslevel: Optional[int] = 9, indent: int = 0):
    """
    Export a dictionary into a (compressed) JSON file.
    :param data: Dictionary with extracted sORFs
    :param file_path: Path to the output file
    :param single_line: write the data into one line
    :param threads: Number of threads for file compression
    :param compresslevel: Compression level for file compression
    :param indent: Indentation of JSON file
    """
    with xopen(file_path, mode='w', compresslevel=compresslevel, threads=threads) as f:
        if single_line:
            json.dump(data, f)
        else:
            f.write(json.dumps(data, indent=indent))


def dict_to_fasta(data: dict[str, str], file_path: Union[str, Path], threads: Optional[int] = None,
                  compresslevel: int = 9):
    """
    Export a dictionary with header: sequence entries to a FASTA file.
    :param data: Dictionary with header: sequence entries
    :param file_path: Path to the output file
    :param threads: Number of threads for file compression
    :param compresslevel: Compression level for file compression
    """
    with xopen(file_path, mode='w', compresslevel=compresslevel, threads=threads) as f:
        for header, sequence in data.items():
            f.write(f">{header}\n{sequence}\n")


def list_to_fasta(data: list[str], file_path: Union[str, Path], mode: str = 'sequences', threads: Optional[int] = None,
                  compresslevel: int = 9):
    """
    Export a dictionary with header: sequence entries to a FASTA file.
    :param data: List with either header and sequences or sequences only
    :param file_path: Path to the output file
    :param mode: 'sequences' or 'interleaved' data contains either only sequences or header and sequences
    :param threads: Number of threads for file compression
    :param compresslevel: Compression level for file compression
    """
    with xopen(file_path, mode='w', compresslevel=compresslevel, threads=threads) as f:
        if mode == 'sequences':
            for header, sequence in enumerate(data):
                f.write(f">{header}\n{sequence}\n")
        elif mode == 'interleaved':
            for i, element in enumerate(data):
                if i == 0 or i % 2 == 0:
                    f.write(f">{element}\n")  # header
                else:
                    f.write(f"{element}\n")  # sequence


def parse_fasta(file_path: Path) -> Generator[tuple[str, str], None, None]:
    """
    Parse a fasta file and yield the headers and sequences.
    :param file_path:
    :return: (header, sequence)
    """
    header: Optional[str] = None
    seqs: list[str] = []
    with xopen(file_path, mode='r') as f:
        for line in f:
            line: str = line.strip()
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seqs)
                    seqs = []
                header = line.lstrip('>')
            else:
                seqs.append(line)
        if header is not None and len(seqs) > 0:
            yield header, ''.join(seqs)


def find_files(path: Union[str, Path], pattern: str, regex: bool = True) -> Generator[Path, None, None]:
    """
    Find all files matching the given pattern in the given directory.
    :param path: Path to directory
    :param pattern: Pattern for file search
    :param regex: Use pattern as regex; else search for substring in filenames
    :return: List of found files
    """
    path = str(Path(path))
    if regex:
        import re
        pattern = re.compile(pattern)

    for root, directories, files in os.walk(path):
        for file in files:
            if regex:
                if pattern.match(file) is not None:
                    yield Path(os.path.join(root, file))
            elif pattern in str(file):
                yield Path(os.path.join(root, file))


# EOF
