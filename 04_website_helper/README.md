# sORFdb - Web helper

## About
This directory contains helper scripts to prepare the upload of the sORFdb database to the server.

## Usage:

Add the cluster information to the database tar archive.
```commandline
add_cluster_to_tar_db.py -c /path/to/db/sorfdb/sorfdb.clusters.1.0.tsv.gz -a /path/to/db/sorfdb/sorfdb.1.0.tar.gz -o /path/to/db/sorfdb/sorfdb.1.0.withClusters.tar.gz 
```

Compute the statistics for small protein families and store them in an additional tar archive for the Elasticsearch cluster.
```commandline
compute_cluster_statistics_for_website.py -c /path/to/db/sorfdb/sorfdb.clusters.1.0.tsv.gz -a /path/to/db/sorfdb/sorfdb.1.0.withClusters.tar.gz -o /path/to/db/sorfdb/sorfdb.1.0.clusterStatistic.tar.gz -m /path/to/db/sorfdb/sorfdb.1.0.hmm.gz -l /path/to/db/sorfdb/raw/alignments 
```
