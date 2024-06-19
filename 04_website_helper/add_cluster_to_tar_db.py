#!/usr/bin/env python3

import gc
import csv
import json
import tarfile
from io import BytesIO
from pathlib import Path
from argparse import ArgumentParser

from xopen import xopen


def read_clusters(cluster_file: Path) -> dict:
    clusters: dict = dict()
    with xopen(cluster_file, 'r') as f:
        tsv_reader = csv.reader(f, delimiter='\t')
        for row in tsv_reader:
            clusters.setdefault(row[0], row[1])
    return clusters


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Export all sORFs to JSON datapoints for Elasticsearch and to TSV.')
    parser.add_argument('--archive', '-a', type=Path,
                        help='Input sORFdb as compressed tar archive')
    parser.add_argument('--clusters', '-c', type=Path,
                        help='Input sORFdb small protein clusters as compressed TSV')
    parser.add_argument('--output', '-o', type=Path, default=Path('./'),
                        help='Output sORFdb tar archive with added cluster information.')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    print('Parse clusters...')
    clusters: dict = read_clusters(args.clusters)

    print('Update sORFdb tar archive...')
    with tarfile.open(args.archive, "r:gz") as tar_in, tarfile.open(args.output, 'w:gz') as tar_out:
        for i, entry in enumerate(tar_in):
            if entry.isfile():
                fh = tar_in.extractfile(entry)
                datapoint: dict = json.load(fh)
                if datapoint['protein'] in clusters:
                    datapoint.setdefault('cluster', clusters[datapoint['protein']])
                    datapoint_string = json.dumps(datapoint, indent=0)
                    del fh, datapoint
                    info = tarfile.TarInfo(name=entry.name)
                    info.mtime = entry.mtime
                    info.mode = entry.mode
                    info.size = len(datapoint_string)
                    tar_out.addfile(tarinfo=info, fileobj=BytesIO(bytes(datapoint_string, 'utf-8')))
            if i % 100000 == 0:
                gc.collect()
                print(f"Processed {i}", flush=True)


if __name__ == '__main__':
    main()

# EOF
