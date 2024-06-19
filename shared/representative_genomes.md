# Instructions

Download representative genomes from:

https://www.ncbi.nlm.nih.gov/assembly?term=%22Bacteria%22%5BOrganism%5D%20AND%20%22representative%20genome%22%5Brefseq%20category%5D%20AND%20%28bacteria%5Bfilter%5D%20AND%20%22latest%20genbank%22%5Bfilter%5D%20AND%20%28%22complete%20genome%22%5Bfilter%5D%20OR%20%22chromosome%20level%22%5Bfilter%5D%20OR%20%22scaffold%20level%22%5Bfilter%5D%29%20AND%20all%5Bfilter%5D%20NOT%20anomalous%5Bfilter%5D%29%20AND%20%28%22latest%20genbank%22%5Bfilter%5D%20AND%20all%5Bfilter%5D%20NOT%20anomalous%5Bfilter%5D%29&cmd=DetailsSearch
https://www.ncbi.nlm.nih.gov/assembly?term=%22Bacteria%22%5BOrganism%5D%20AND%20%22representative%20genome%22%5Brefseq%20category%5D%20AND%20%28bacteria%5Bfilter%5D%20AND%20%22latest%20genbank%22%5Bfilter%5D%20AND%20%28%22complete%20genome%22%5Bfilter%5D%20OR%20%22chromosome%20level%22%5Bfilter%5D%20OR%20%22scaffold%20level%22%5Bfilter%5D%29%20AND%20all%5Bfilter%5D%20NOT%20anomalous%5Bfilter%5D%29&cmd=DetailsSearch

Use "Download Assemblies" -> Source database: GenBank -> File type: Assembly statistics report (.txt)

Untar

Run `ls GCA_*.txt | cut -d '.' -f 1 | sort | uniq > representative_genomes.txt`
