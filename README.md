![x86-64](https://img.shields.io/static/v1?label=%E2%80%8B&message=x86-64&color=yellow&logo=PCGamingWiki&logoColor=white)
![Linux](https://img.shields.io/static/v1?label=%E2%80%8B&message=Linux&color=00A98F&logo=linux&logoColor=white)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-brightgreen.svg)](https://github.com/ag-computational-bio/sorfdb/blob/main/LICENSE)


# sORFdb - A database for sORFs, small proteins, and small protein families in bacteria

sORFdb is a comprehensive, taxonomically independent database dedicated to **sORF** and **small protein** sequences, 
**small protein families** and related information in bacteria. It aims to improve the findability and classification of 
sORFs, small proteins, and their functions in bacteria, thereby supporting their future detection and consistent 
annotation.

The website of sORFdb is available at https://sorfdb.computational.bio/.

The database can be dowloaded from [Zenodo](https://doi.org/10.5281/zenodo.10688271)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10688271.svg)](https://doi.org/10.5281/zenodo.10688271)


## Contents

- [Description](#description)
- [Installation](#installation)

## Description

This repository contains the workflows used to create the sORFdb database starting from the data aggregation 
(`01_download`), the data processing (`02_processing`), the clustering of small proteins and identification of small 
protein families (`03_autoclust`) and helper scripts to prepare the database for the server (`04_website-helper`). The 
Jupyter notebook used for conducting the analysis of the data for the manuscript is also available (`05_analysis`).

To create the database from scratch, it is highly recommended to have access to a SLURM cluster or expand the nextflow 
config files with another nextflow executor.

## Installation

Clone this GitHub repository to your local system. Run the according Nextflow scripts in the subdirectories. Their 
usage is described in their respective README files.

**Requirements**

- [Nextflow](https://nextflow.io/) 
- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/en/latest/index.html)  


