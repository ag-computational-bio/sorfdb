#!/usr/bin/env python3

import os
from pathlib import Path
from argparse import ArgumentParser
from typing import Union

import polars as pl


def parse_clm(file_path: Path) -> pl.DataFrame:
    """
    Read the clm info file content and convert it into a Dataframe.
    :param file_path: path to clminfo.txt file
    :return: Dataframe with clm info content
    """
    clm_data: dict[str, list[Union[str, float]]] = {'inflation': []}

    float_keys: set[str] = {'eff', 'mod', 'mf', 'af', 'ctr', 'avg'}
    int_keys: set[str] = {'ncl', 'max', 'min', 'DGI', 'TWI', 'TWL', 'sgl', 'qrt'}

    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('eff'):
                fields: list[str] = [l.strip() for l in line.strip().split() if len(l) > 0]
                for field in fields:
                    key, value = field.split('=')
                    clm_data.setdefault(key, [])
                    if key in float_keys:
                        clm_data[key].append(float(value))
                    elif key in int_keys:
                        clm_data[key].append(int(value))
                    else:
                        clm_data[key].append(value)
                    if key == 'src':
                        tmp_inflation: str = value.split('.')[-1].lstrip('I')
                        inflation: float = float(tmp_inflation[0] + '.' + tmp_inflation[1:])
                        clm_data['inflation'].append(inflation)

    clm_df: pl.DataFrame = pl.from_dict(clm_data)
    return clm_df.sort('inflation')


def max_modularity(clm_df: pl.DataFrame) -> pl.DataFrame:
    inflation_df: pl.DataFrame = clm_df.filter(
        pl.col('mod') == pl.col('mod').max()
    ).filter(
        pl.col('eff') == pl.col('eff').max()
    ).filter(
        pl.col('inflation') == pl.col('inflation').min()
    )
    return inflation_df


def max_efficacy(clm_df: pl.DataFrame) -> pl.DataFrame:
    inflation_df: pl.DataFrame = clm_df.filter(
        pl.col('eff') == pl.col('eff').max()
    ).filter(
        pl.col('mod') == pl.col('mod').max()
    ).filter(
        pl.col('inflation') == pl.col('inflation').min()
    )
    return inflation_df


def max_combined(clm_df: pl.DataFrame) -> pl.DataFrame:
    inflation_modularity_df: pl.DataFrame = max_modularity(clm_df)
    inflation_efficacy_df: pl.DataFrame = max_efficacy(clm_df)
    inflation_modularity: float = inflation_modularity_df.select(pl.col('inflation')).to_series().to_list()[0]
    inflation_efficacy: float = inflation_efficacy_df.select(pl.col('inflation')).to_series().to_list()[0]
    inflation: float = round((inflation_modularity + inflation_efficacy) / 2, 1)
    inflation_df: pl.DataFrame = clm_df.filter(
        pl.col('inflation') == inflation
    )
    return inflation_df


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Import a mcl abc file and tab files with clusters and export them for use in cytoscape.')
    parser.add_argument('--input', '-i', type=Path, help='Input clm file for MCL clustering')
    parser.add_argument('--id', '-d', type=str, help='Id of the subgraph.')
    parser.add_argument('--metric', '-m', type=str, default='mod',
                        help='Metric to select the best MCL inflation parameter. Use Efficacy (eff) or Modularity (mod). (default=mod)')
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()
    pl.Config.set_tbl_rows(100)
    pl.Config.set_tbl_cols(16)
    clm_df: pl.DataFrame = parse_clm(args.input)

    if args.metric == 'mod':
        inflation_df: pl.DataFrame = max_modularity(clm_df)
    elif args.metric == 'eff':
        inflation_df: pl.DataFrame = max_efficacy(clm_df)
    elif args.metric == 'comb':
        inflation_df: pl.DataFrame = max_combined(clm_df)
    else:
        print(f'Metric {args.metric} is not supported. Use Efficacy (eff) or Modularity (mod).')
        exit(1)

    inflation: float = inflation_df.select(pl.col('inflation')).to_series().to_list()[0]
    print(inflation)
    inflation_file: str = inflation_df.select(pl.col('src')).to_series().to_list()[0]
    base_path: Path = Path('/'.join(str(args.input).split('/')[:-1]))
    os.symlink(base_path.joinpath(inflation_file), base_path.joinpath(f'{args.id}.mci'))


if __name__ == '__main__':
    main()

# EOF
