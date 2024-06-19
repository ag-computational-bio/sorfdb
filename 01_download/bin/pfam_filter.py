#!/usr/bin/env python3

import re
import sys
from pathlib import Path

from xopen import xopen


def profile_splitter(pfam_path: Path, out_path: Path):
    with xopen(pfam_path, 'r') as f:
        file_contents = f.read().split('//')[0:-1]

    pattern = re.compile(r'\nLENG\s+([0-9]|[1-9][0-9]|100|101)\n')

    out_file: Path = out_path.joinpath('pfam_short')

    with open(out_file, "w") as f:
        for profile in file_contents:
            match = re.search(pattern, profile)
            if match is not None:
                if profile[:1] == '\n':
                    profile = profile[1:]
                f.write(profile + "//\n")


def main():
    _, pfam_path, out_path = sys.argv
    pfam_path: Path = Path(pfam_path)
    out_path: Path = Path(out_path)

    profile_splitter(pfam_path, out_path)


if __name__ == '__main__':
    main()

# EOF
