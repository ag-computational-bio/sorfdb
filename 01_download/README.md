# sORFdb - Data aggregation workflow

## Description

This workflow collects all the data needed for the sORFdb database. All data for the processing steps are stored 
in their respective subdirectories in the directory specified by `--db`. 

## Usage:
```commandline
nextflow run db-setup.nf --db /path/to/db/ --workDir /working_directory/ --mamba true
```