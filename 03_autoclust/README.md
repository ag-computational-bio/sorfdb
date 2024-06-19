# sORFdb - Autoclust

## About
Autoclust is a graph based clustering approach using different heuristics and pruning strategies before using MCL for 
the clustering. It uses normalized bit scores in the form of score ratio values (SRV) as edge weights.

## Installation
Autoclust is designed for Linux systems and can be run locally or on a SLURM cluster.

## Usage:
```commandline
nextflow run -profile slurm autoclust.nf --db /path/to/db/ --mamba true --input /path/to/db/sorfdb/sprot.clustering.faa.gz --cov 70
```