#!/usr/bin/env python3

from math import floor
from pathlib import Path
from argparse import ArgumentParser

import numpy as np
import polars as pl


def separate_graphs(graph_df: pl.DataFrame) -> list[pl.DataFrame]:
    """
    Separate completely unlinked graphs.
    :param graph_df:
    :return: list of unlinked subgraphs
    """
    subgraphs: list[pl.DataFrame] = []
    subgraph_nodes: pl.Series
    subgraph: pl.DataFrame = pl.from_dict({'query': [], 'subject': [], 'srv': []})
    while graph_df.shape[0] > 0:
        out_node: str = graph_df.select(pl.col('query')).head(n=1).to_series().to_list()[0]
        linked_nodes: set[str] = set(graph_df.filter(
            pl.col('query') == out_node
        ).select(pl.col('subject')).to_series().to_list())
        linked_nodes.add(out_node)

        node_count: int = 0
        while len(linked_nodes) > node_count:
            node_count = len(linked_nodes)
            subgraph_nodes = pl.Series('subnodes', tuple(linked_nodes))
            subgraph = graph_df.filter(
                (pl.col('query').is_in(subgraph_nodes)) |
                (pl.col('subject').is_in(subgraph_nodes))
            )
            linked_nodes = set(subgraph.select(pl.col('query')).unique().to_series().to_list() +
                               subgraph.select(pl.col('subject')).unique().to_series().to_list())

        subgraphs.append(
            subgraph.rechunk()
        )
        subgraph_nodes = pl.Series('subnodes', tuple(linked_nodes))
        graph_df = graph_df.filter(
            ~(pl.col('query').is_in(subgraph_nodes)) &
            ~(pl.col('subject').is_in(subgraph_nodes))
        ).rechunk()
    return subgraphs


def mean_node_degree(graph: pl.DataFrame) -> float:
    """
    Calculate the mean node degree within a subgraph.
    :param graph: graph
    :return: mean node degree
    """
    degree: float = graph.select(
        pl.col('query')
    ).groupby('query').count().select(
        pl.col('count')
    ).mean().to_series().to_list()[0]
    return degree


def calculate_bins(node_degrees: np.ndarray) -> np.ndarray:
    """
    Calculate bin sizes on the mean node degree of the subgraphs.
    Choose the option with the most bins.
    :param node_degrees: array with mean node degrees
    :return: array with bin ranges
    """
    bins_doane: np.ndarray = np.histogram_bin_edges(node_degrees, bins='doane')
    bins_fd: np.ndarray = np.histogram_bin_edges(node_degrees, bins='fd')
    if len(bins_doane) > len(bins_fd):
        return bins_doane
    return bins_fd


def heuristic_threshold(subgraphs: list[pl.DataFrame], min_srv: int) -> list[pl.DataFrame]:
    """
    Apply the heuristical threshold proposed by Leonard Apeltsin et al., 2011, https://doi.org/10.1093/bioinformatics/btq655.
    Return subgraphs with applied computation.
    :param subgraphs:
    :param min_srv:
    :return:
    """
    filtered_subgraphs: list[pl.DataFrame] = []
    for subgraph in subgraphs:
        if subgraph.shape[0] <= 3:  # new threshold would produce singletons
            filtered_subgraphs.append(subgraph)
            continue
        min_edge_weight: int = int(floor(subgraph.select(pl.col('srv')).min().to_series().to_list()[0]))
        min_th: int = max(min_srv, min_edge_weight)
        for th in range(min_th, 100):
            # Nn(Th), is the number of nodes connected by one or more edges at threshold Th
            # SE(Th), the number edges remaining after threshold Th is applied
            # Nsv(Th) = SE(Th)/Nn(Th); equivalent to the average weighted node degree at threshold Th
            # threshold estimation heuristic: b is approximate to the minimum threshold Th at which dNsv(Th)/dTh > 0.
            # If Nsv does not increase at any point in the distribution, then no threshold is returned.
            graph_th: pl.DataFrame = subgraph.filter(pl.col('srv') >= th)
            nn: int = graph_th.select(pl.col('query')).vstack(graph_th.select(pl.col('subject').alias('query'))).unique().shape[0]
            se: int = graph_th.shape[0]
            try:
                nsv: float = se / nn
            except ZeroDivisionError:
                print(f'Could not apply new threshold={th}.')
                print(subgraph)
                nsv = 0.0

            if nsv > 1.0:
                if th > min_edge_weight:  # new higher th
                    filtered_subgraphs.extend(separate_graphs(graph_th))
                    break
                else:  # Nsv(Th) is already > 0
                    filtered_subgraphs.append(subgraph)
                    break
        else:  # no Nsv(Th) > 0
            filtered_subgraphs.append(subgraph)

    return filtered_subgraphs


def combine_subgraphs(subgraphs: list[pl.DataFrame]) -> list[pl.DataFrame]:
    """
    Combine subgraphs with similar mean node degrees.
    :param subgraphs:
    :return:
    """
    mean_node_degrees: np.ndarray = np.array([mean_node_degree(graph) for graph in subgraphs])
    mean_edge_weigths: np.ndarray = np.array([graph.select(pl.col('srv')).mean().to_series().to_list()[0] for graph in subgraphs])
    graph_criterion: np.ndarray = mean_node_degrees * mean_edge_weigths
    bins: np.ndarray = calculate_bins(graph_criterion)
    assigned_bins: np.ndarray = np.digitize(graph_criterion, bins)

    binned_subgraphs: list[list[pl.DataFrame]] = [[] for i in range(0, len(bins))]
    for i, subgraph in enumerate(subgraphs):
        graph_bin: int = assigned_bins[i] - 1  # has one based index
        binned_subgraphs[graph_bin].append(subgraph)

    combined_subgraphs: list[pl.DataFrame] = []
    for binned_graphs in binned_subgraphs:
        if len(binned_graphs) == 0:
            continue
        combined_subgraph: pl.DataFrame = binned_graphs[0]
        if len(binned_graphs) > 1:
            for subgraph in binned_graphs[1:]:
                combined_subgraph.vstack(subgraph, in_place=True)
        combined_subgraphs.append(
            combined_subgraph
        )
    del binned_subgraphs
    return [graph.rechunk() for graph in combined_subgraphs]


def parse_arguments():
    """
    Argument parser
    :return: Command line arguments.
    """
    parser = ArgumentParser(description='Import a mcl abc file and tab files with clusters and export them for use in cytoscape.')
    parser.add_argument('--abc', '-a', type=Path, help='Input mcl tab file with clusters')
    parser.add_argument('--threshold', '-d', type=int, default=1, help='SRV minimum threshold')
    parser.add_argument('--output', '-o', type=Path, default=Path('./'), help='Output folder path')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    print('Import graph...')
    grapf_df: pl.DataFrame = pl.scan_csv(
        args.abc,
        separator='\t',
        has_header=False,
        new_columns=['query', 'subject', 'srv']
    ).collect()

    print('Separate graphs...')
    subgraphs: list[pl.DataFrame] = separate_graphs(grapf_df)
    del grapf_df

    print('Apply heuristic...')
    subgraphs = heuristic_threshold(subgraphs, 30)

    print('Bin graphs...')
    subgraphs = combine_subgraphs(subgraphs)

    print('Write subgraphs...')
    for i, subgraph in enumerate(subgraphs):
        subgraph.write_csv(args.output.joinpath(f'{i}.subgraph.tsv'),
                           separator='\t',
                           has_header=False)


if __name__ == '__main__':
    main()

# EOF
