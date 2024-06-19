# sORFdb - Data processing workflow

## Description

This workflow processes all collected data stored in the directory specified by `--db` for the sORFdb database. The 
files needed for the webserver are stored in `--db /dbPath/sorfdb/`.

If you want to use the workflow on your system, change the SLURM executor in the `nextflow.config` according to your 
SLURM setup.

## Usage:
```commandline
nextflow run -profile slurm db-processing.nf --db /path/to/db/ --workDir /working_directory/ --mamba true
```