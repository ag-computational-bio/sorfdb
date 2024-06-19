#!/usr/bin/env python3

import re

DNA_PATTERN: re.Pattern = re.compile(r'[ACGT]+', re.IGNORECASE)
DNA_AMBIGUOUS_PATTERN: re.Pattern = re.compile(r'[ACGTUWSMKRYBDHVN]+', re.IGNORECASE)
PROTEIN_PATTERN: re.Pattern = re.compile(r'[ACDEFGHIKLMNPQRSTVWY]+', re.IGNORECASE)
PROTEIN_EXTENDED_PATTERN: re.Pattern = re.compile(r'[ACDEFGHIKLMNOPQRSTUVWY]+', re.IGNORECASE)
PROTEIN_AMBIGUOUS_PATTERN: re.Pattern = re.compile(r'[ACDEFGHIKLMNPQRSTVWXY]+', re.IGNORECASE)
PROTEIN_EXTENDED_AMBIGUOUS_PATTERN: re.Pattern = re.compile(r'[ACDEFGHIKLMNOPQRSTUVWXY]+', re.IGNORECASE)


def is_dna(seq: str, pattern: re.Pattern = DNA_PATTERN) -> bool:
    """
    Check if the DNA sequence contains only nucleotides from a given pattern.
    :param seq: DNA sequence
    :param pattern: DNA patter to match (default ACGT)
    :return: Boolean
    """
    if pattern.fullmatch(seq) is None:
        return False
    return True


def is_protein(seq: str, pattern: re.Pattern = PROTEIN_PATTERN) -> bool:
    """
    Check if the protein sequence contains only amino acids from a given pattern.
    :param seq: protein sequence
    :param pattern: protein patter to match (default unambiguous)
    :return: Boolean
    """
    if pattern.fullmatch(seq) is None:
        return False
    return True

# EOF
