#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import static groovy.io.FileType.*
import static groovy.io.FileVisitResult.*

/*
 * If the 'output' is not empty print a warning.
 * You might override files if the 'output' is not empty.
 */

output = file(params.db)

if (output != null && output.list().size() != 0) {
    println("Warning: Directory $output is not empty. You may override existing files.")
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////WORKFLOW///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow {
    query_Genbank()
    download_Genbank(query_Genbank.out.genomes.flatten())
    download_UniRef100()
    extract_UniRef100(download_UniRef100.out.uniref, download_UniRef100.out.taxonomy)
    download_SwissProt()
    extract_SwissProt(download_SwissProt.out.complete, extract_UniRef100.out.uniref)
    download_UniProt_filtered()
    extract_UniProt_filtered(download_UniProt_filtered.out.complete, extract_UniRef100.out.uniref)
    download_AntiFam()
    download_Pfam()
    setup_expert_DBs(extract_SwissProt.out.faa,
                     extract_UniProt_filtered.out.faa)
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////GENBANK////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Query for Genbank database files.
 */

process query_Genbank {
    publishDir "${params.db}/log", mode: 'copy', pattern: "*.log"

    conda "$baseDir/envs/db-setup.yml"

    when:
        params.update == true

    output:
        path("genomes_*"), emit: genomes
        path("*.log"), emit: log

    script:
        """
        wget https://ftp.ncbi.nih.gov/genbank/GB_Release_Number -O GB_Release_Number.log

        curl 'https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt' | \
        gawk -F "\t" '(\$12=="Complete Genome" || \$12=="Chromosome" || \$12=="Scaffold") && \$11=="latest"{print \$20}' | \
        gawk 'BEGIN{FS=OFS="/";filesuffix="genomic.gbff.gz"}{ftpdir=\$0;asm=\$10;file=asm"_"filesuffix;print ftpdir,file}' | \
        gawk '{if (\$0 != "na/_genomic.gbff.gz") print \$0}' | \
        gawk -F ':' '{print "https:"\$2}' > all_genomes

        date > timestamp.log
        cat all_genomes | wc -l > genbank_genomes.log
        CHUNKSIZE=\$(cat genbank_genomes.log | gawk '{print int(\$1/40)}')

        split -l \$CHUNKSIZE all_genomes genomes_
        """
}

/*
 * Download and set up the Genbank database.
 */

process download_Genbank {
    publishDir "${params.db}/genbank", mode: 'move', pattern: "*genomic.gbff.gz"

    errorStrategy { sleep(Math.pow(2, task.attempt) * 100000 as long); return 'retry' }
    maxRetries = 5
    maxForks (params.threads >= 4 ? 4 : params.threads)
    conda "$baseDir/envs/db-setup.yml"

    input:
        file(chunk)

    output:
        path("*genomic.gbff.gz")

    script:
        """
        wget --quiet --continue -i ${chunk}
        """
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////Swiss/UniProt/////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
 * Download the UniRef100 database.
 */

process download_UniRef100 {
    publishDir "${params.db}/proteinDB", mode: 'copy', pattern: "bacteria_tax_ids.tsv.gz"
    publishDir "${params.db}/log", mode: 'copy', pattern: "relnotes.txt"
    publishDir "${params.db}/log", mode: 'copy', pattern: "*.log"
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries = 3
    conda "$baseDir/envs/db-setup.yml"

    output:
        path("uniref100.xml.gz"), emit: uniref
        path("bacteria_tax_ids.tsv.gz"), emit: taxonomy
        path("relnotes.txt"), emit: relnotes

    script:
        """
        wget https://ftp.expasy.org/databases/uniprot/current_release/relnotes.txt &

        curl --connect-timeout 30 --retry 3 --retry-delay 5 --globoff 'https://rest.uniprot.org/taxonomy/stream?compressed=true&fields=id%2Cscientific_name%2Clineage&format=tsv&query=%28ancestor%3A2%29' \
        > bacteria_tax_ids.tsv.gz &

        wget https://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref100/uniref100.xml.gz
        """
}

/*
 * Extract all bacterial small protein entries from the UniRef100 database.
 */

process extract_UniRef100 {
    publishDir "${params.db}/log", mode: 'copy', pattern: "*.log"
    memory { 32.GB * task.attempt * Math.pow(2, task.attempt) }
    maxRetries = 3
    conda "$baseDir/envs/db-setup.yml"

    input:
        path(uniref100)
        path(taxonomy)

    output:
        path("uniref100.faa.gz"), emit: uniref
        path("*.log")

    script:
        """
        uniref100_parser.py --uniref $uniref100 --taxonomy $taxonomy \
        --output ./uniref100.faa.gz --threads ${task.cpus}

        ln -s uniref100.faa.gz uniprot.all.faa.gz
        zgrep -c '>' uniprot.all.faa.gz > uniprot.all.log
        """
}

/*
 * Download the Swiss-Prot database.
 */

process download_SwissProt {
    publishDir "${params.db}/proteinDB", mode: 'copy', pattern: "*.faa.gz"
    publishDir "${params.db}/log", mode: 'copy', pattern: "*.log"

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries = 3
    conda "$baseDir/envs/db-setup.yml"

    output:
        path("swissprot.complete.faa.gz"), emit: complete
        path("swissprot.complete.log"), emit: log

    script:
        """
        curl --connect-timeout 30 --retry 3 --retry-delay 5 --globoff \
        'https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28taxonomy_id%3A2%29%20AND%20%28length%3A%5B%2A%20TO%20${params.length}%5D%29%20AND%20%28fragment%3Afalse%29%29%20AND%20%28reviewed%3Atrue%29' \
        > swissprot.complete.faa.gz

        zgrep -c '>' swissprot.complete.faa.gz > swissprot.complete.log
        """
}

/*
 * Extract the UniRef100 entries of the Swiss-Sprot database entries.
 */

process extract_SwissProt {
    publishDir "${params.db}/proteinDB", mode: 'copy', pattern: "*.faa.gz"
    publishDir "${params.db}/log", mode: 'copy', pattern: "*.log"

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries = 3
    conda "$baseDir/envs/db-setup.yml"

    input:
        path(swissprot)
        path(uniref100)

    output:
        path("swissprot.faa.gz"), emit: faa
        path("swissprot.log"), emit: log

    script:
        """
        extract_uniref_entries.py --input $swissprot --uniref $uniref100 --prefix swissprot \
        --output swissprot.faa.gz --threads ${task.cpus}

        zgrep -c '>' swissprot.faa.gz > swissprot.log
        """
}

/*
 * Download the UniProt database filtered for evidence.
 */

process download_UniProt_filtered {
    publishDir "${params.db}/proteinDB", mode: 'copy', pattern: "*.faa.gz"
    publishDir "${params.db}/log", mode: 'copy', pattern: "*.log"

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries = 3
    conda "$baseDir/envs/db-setup.yml"

    output:
        path("uniprot.filtered.complete.faa.gz"), emit: complete
        path("*.log"), emit: log

    script:
        """
        curl --connect-timeout 30 --retry 3 --retry-delay 5 --globoff \
        'https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28taxonomy_id%3A2%29%20AND%20%28length%3A%5B%2A%20TO%20${params.length}%5D%29%20AND%20%28fragment%3Afalse%29%20AND%20%28%28existence%3A1%29%20OR%20%28existence%3A2%29%20OR%20%28existence%3A3%29%29%29' \
        > uniprot.filtered.complete.faa.gz

        zgrep -c '>' uniprot.filtered.complete.faa.gz > uniprot.filtered.complete.log
        """
}

/*
 * Extract the UniRef100 entries of the evidence filtered database entries.
 */

process extract_UniProt_filtered {
    publishDir "${params.db}/proteinDB", mode: 'copy', pattern: "*.faa.gz"
    publishDir "${params.db}/log", mode: 'copy', pattern: "*.log"

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries = 3
    conda "$baseDir/envs/db-setup.yml"

    input:
        path(uniprot)
        path(uniref100)

    output:
        path("uniprot.filtered.faa.gz"), emit: faa
        path("*.log"), emit: log

    script:
        """
        extract_uniref_entries.py --input $uniprot --uniref $uniref100 --prefix uniprot \
        --output uniprot.filtered.faa.gz --threads ${task.cpus}

        zgrep -c '>' uniprot.filtered.faa.gz > uniprot.filtered.log
        """
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////AntiFam////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
 * Download the current AntiFam database.
 */

process download_AntiFam {
    publishDir "${params.db}/antifam", mode: 'copy', pattern: "antifam.h3*"

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries = 3
    conda "$baseDir/envs/db-setup.yml"

    output:
        file("antifam.h3*")

    script:
        """
        mkdir antifam-dir
        cd antifam-dir
        curl -O "https://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz"
        tar -xzf Antifam.tar.gz
        cd ..
        mv antifam-dir/AntiFam_Bacteria.hmm antifam
        hmmpress antifam
        rm -r antifam antifam-dir/
        """
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////Pfam////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
 * Download the current Pfam database.
 */

process download_Pfam {
    publishDir "${params.db}/pfam", mode: 'copy', pattern: "pfam_short.h3*"

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries = 3
    conda "$baseDir/envs/db-setup.yml"

    output:
        file("pfam_short.h3*")

    script:
        """
        curl -O "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
        pfam_filter.py Pfam-A.hmm.gz ./
        hmmpress pfam_short
        rm Pfam-A.hmm* pfam_short
        """
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////smProt////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Create a diamond database from small proteins from all databases.
 */

process setup_expert_DBs {
    publishDir "${params.db}/proteinDB", mode: 'copy', pattern: "*.faa.gz"
    publishDir "${params.db}/proteinDB", mode: 'copy', pattern: "*.dmnd"
    publishDir "${params.db}/log", mode: 'copy', pattern: "*.log"

    conda "$baseDir/envs/db-setup.yml"

    input:
        file("swissprot.faa.gz")
        file("uniprot.filtered.faa.gz")

    output:
        file("*.faa.gz")
        file("*.dmnd")
        file("*.log")

    script:
        """
        # Order: low to high
        # FILTERED
        unique_faa_merger.py filtered ${baseDir}/data/${params.length}/smprot_ribo_seq.faa.gz \
        ${baseDir}/data/${params.length}/smprot_kdb.faa.gz ${baseDir}/data/${params.length}/smprot_literature.faa.gz \
        uniprot.filtered.faa.gz swissprot.faa.gz

        zgrep -c '>' combined_sorf_db.filtered.faa.gz > sORFs.filtered.log

        diamond makedb --in combined_sorf_db.filtered.faa.gz --db sORFs.filtered.dmnd

        gzip -dck combined_sorf_db.filtered.faa.gz | bioawk -c fastx '{if(length(\$seq)<=50){print ">" \$name ORS \$seq}}' | \
        gzip -c > combined_sorf_db.filtered.50.faa.gz
        diamond makedb --in combined_sorf_db.filtered.50.faa.gz --db sORFs.filtered.50.dmnd
        """
}

// EOF
