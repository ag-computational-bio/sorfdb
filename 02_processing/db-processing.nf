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

/*
 * Execution information
 */

println("Executor=" + params.executor)
println("Threads=" + params.threads)
params.mamba ? println("Using Mamba.") : println("Using Conda.")
println("\n")


///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////WORKFLOW///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow {
    ch_input_genbank = Channel.fromPath("${params.db}/genbank/*genomic.gbff.gz")
                                .map{x -> {
                                     def blocks = x.name =~ /(.*)_genomic\.gbff\.gz/
                                     return [blocks[0][1].split('_')[0..1].join('_'), x]
                                }}
    
    parse_Genbank_2_json(ch_input_genbank)

    int annotated_counter = 0
    precollect_annotated(parse_Genbank_2_json.out.annotated.toSortedList().flatten().buffer(size: 1000,
                                                                                            remainder: true).map{[annotated_counter++, it]})
    int hypothetical_counter = 0
    precollect_hypotheticals(parse_Genbank_2_json.out.hypothetical.toSortedList().flatten().buffer(size: 1000,
                                                                                                   remainder: true).map{[hypothetical_counter++, it]})

    int blastx_counter = 0
    precollect_blastx(parse_Genbank_2_json.out.blastx.toSortedList().flatten().buffer(size: 1000,
                                                                                      remainder: true).map{[blastx_counter++, it]})

    concat_annotated(precollect_annotated.out.collect())
    concat_hypotheticals(precollect_hypotheticals.out.collect())
    concat_blastx(precollect_blastx.out.collect())

    hmmer_annotated(concat_annotated.out)
    hmmer_hypotheticals(concat_hypotheticals.out)
    hmmer_blastx(concat_blastx.out)

    verify_hypotheticals(hmmer_hypotheticals.out)
    verify_blastx(hmmer_blastx.out)

    deduplicate_valid(hmmer_annotated.out.valid, verify_hypotheticals.out)

    ch_input_swissprot = Channel.fromPath("${params.db}/proteinDB/swissprot.complete.faa.gz")
    ch_input_uniprot = Channel.fromPath("${params.db}/proteinDB/uniprot.filtered.complete.faa.gz")
    ch_input_smprot = Channel.fromPath("${baseDir}/data/${params.length}/smProt_combined_database.tsv")
    ch_input_taxonomy = Channel.fromPath("${params.db}/proteinDB/bacteria_tax_ids.tsv.gz")
    sorfdb(deduplicate_valid.out.all, verify_blastx.out, ch_input_swissprot, ch_input_uniprot, ch_input_smprot, ch_input_taxonomy)
    sorfdb_taxonomy(sorfdb.out)
    blast_export(deduplicate_valid.out.all, verify_blastx.out, ch_input_swissprot, ch_input_uniprot, ch_input_smprot)
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////GENBANK////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Extract all valid and hypothetical as well as possible missing sORFs from genbank annotations.
 */

process parse_Genbank_2_json {
    memory { 512.MB * task.attempt * Math.pow(2, task.attempt) }
    time { 10.min * task.attempt }
    maxRetries = 5

    input:
        tuple val(id), path("${id}.gbff.gz")
    output:
        path("${id}.annotated.json.gz"), optional: true, emit: annotated
        path("${id}.hypothetical.json.gz"), optional: true, emit: hypothetical
        path("${id}.blastx.json.gz"), optional: true, emit: blastx
    script:
        """
        gbff_parser.py --id ${id} --gbff ${id}.gbff.gz --db \$(realpath ${params.db}/proteinDB/sORFs.filtered.50.dmnd) --threads ${task.cpus}
        """
}

/*
 * Pre-collect annotated sORF JSON files.
 */

process precollect_annotated {
    cpus (params.threads >= 4 ? 4 : params.threads)
    memory '2 GB'

    input:
        tuple val(num), path(annotatedJson)
    output:
        path("annotated.${num}.json.gz")
    script:
        """
        pigz -dck *.annotated.json.gz | jq --slurp 'add' | sed 's/^ *//g' | pigz -n -p ${task.cpus} -3 > annotated.${num}.json.gz
        """
}


/*
 * Pre-collect hypothetical JSON files and concatenate them.
 */

process precollect_hypotheticals {
    cpus (params.threads >= 4 ? 4 : params.threads)
    memory '2 GB'

    input:
        tuple val(num), path(hypotheticalJson)
    output:
        path("hypothetical.${num}.json.gz")
    script:
        """
        pigz -dck *.hypothetical.json.gz | jq --slurp 'add' | sed 's/^ *//g' | pigz -n -p ${task.cpus} -3 > hypothetical.${num}.json.gz
        """
}


/*
 * Pre-collect blastx JSON files and concatenate them.
 */

process precollect_blastx {
    cpus (params.threads >= 4 ? 4 : params.threads)
    memory '2 GB'

    input:
        tuple val(num), path(hypotheticalJson)
    output:
        path("blastx.${num}.json.gz")
    script:
        """
        pigz -dck *.blastx.json.gz | jq --slurp 'add' | sed 's/^ *//g' | pigz -n -p ${task.cpus} -3 > blastx.${num}.json.gz
        """
}


/*
 * Concat the annotated sorf entries.
 */

process concat_annotated {
    memory (params.large ? '80 GB' : '24 GB')

    input:
        path(json)

    output:
        path("all.genbank.json")

    script:
        """
        pigz -dck annotated.*.json.gz | jq --slurp 'add' | sed 's/^ *//g' > all.genbank.json
        """
}


/*
 * Concat the hypothetical sorf entries.
 */

process concat_hypotheticals {
    memory (params.large ? '128 GB' : '48 GB')

    input:
        path(json)

    output:
        path("hypothetical.genbank.json")

    script:
        """
        pigz -dck hypothetical.*.json.gz | jq --slurp 'add' | sed 's/^ *//g' > hypothetical.genbank.json
        """
}


/*
 * Concat the blastx sorf entries.
 */

process concat_blastx {
    memory (params.large ? '200 GB' : '48 GB')

    input:
        path(json)

    output:
        path("blastx.json")

    script:
        """
        pigz -dck blastx.*.json.gz | jq --slurp 'add' | sed 's/^ *//g' > blastx.json
        """
}


/*
 * Perform HMMER-searches against antifam and Pfam for annotated sORFs.
 */

process hmmer_annotated {
    cpus params.threads
    memory (params.large ? '300 GB' : '64 GB')

    input:
        path("unfiltered.annotated.json")

    output:
        path("annotated.json"), emit: valid
        path("antifam.annotated.json.gz"), optional: true, emit: antifam

    script:
        """
        hmm.py --data unfiltered.annotated.json --source annotated --antifam ${params.db}/antifam/antifam.h3m \
               --dump_antifam --pfam ${params.db}/pfam/pfam_short.h3m --threads ${task.cpus}
        """
}


/*
 * Perform HMMER-searches against antifam and Pfam for hypothetical sORFs.
 */

process hmmer_hypotheticals {
    queue = params.longqueue
    cpus params.threads
    memory (params.large ? '400 GB' : '64 GB')

    input:
        path("unfiltered.hypothetical.json")

    output:
        tuple path("hypothetical.json"), path("hypothetical.faa.gz")

    script:
        """
        hmm.py --data unfiltered.hypothetical.json --source hypothetical --antifam ${params.db}/antifam/antifam.h3m \
               --remove_antifam --pfam ${params.db}/pfam/pfam_short.h3m --threads ${task.cpus}

        export_faa.py --input hypothetical.json --output ./hypothetical.faa.gz --threads ${task.cpus}
        """
}


/*
 * Perform HMMER-searches against antifam and Pfam for blastx sORFs.
 */

process hmmer_blastx {
    cpus params.threads
    memory (params.large ? '32 GB' : '8 GB')

    input:
        path("unfiltered.blastx.json")

    output:
        tuple path("blastx.json"), path("blastx.faa.gz")

    script:
        """
        hmm.py --data unfiltered.blastx.json --source blastx --antifam ${params.db}/antifam/antifam.h3m \
               --remove_antifam --pfam ${params.db}/pfam/pfam_short.h3m --threads ${task.cpus}

        export_faa.py --input blastx.json --output ./blastx.faa.gz --threads ${task.cpus}
        """
}


/*
 * Verify hypothetical sorfs with homology for known sorfs.
 */

process verify_hypotheticals {
    queue = params.longqueue
    cpus (params.executor != 'local' ? params.maxThreads : params.threads)
    memory (params.large ? '240 GB' : '32 GB')

    input:
        tuple path("hypothetical.genbank.json"), path("hypothetical.genbank.faa.gz")

    output:
        path("verified.genbank.json")

    script:
        """
        validator.py --hypotheticals hypothetical.genbank.json --fasta hypothetical.genbank.faa.gz \
        --db ${params.db}/proteinDB/sORFs.filtered.dmnd --threads ${task.cpus}

        export_faa.py --input verified.genbank.json --output ./hypothetical.faa.gz --threads ${task.cpus}
        """
}


/*
 * Verify blastx sorfs with homology for known sorfs.
 */

process verify_blastx {
    queue = params.longqueue
    cpus params.threads
    memory (params.large ? '32 GB' : '8 GB')

    input:
        tuple path("unverified.blastx.json"), path("unverified.blastx.faa.gz")

    output:
        path("verified.blastx.json")

    script:
        """
        validator.py --hypotheticals unverified.blastx.json --fasta unverified.blastx.faa.gz \
        --db ${params.db}/proteinDB/sORFs.filtered.dmnd --source blastx --threads ${task.cpus}
        """
}


/*
 * Deduplicate all valid sorf entries.
 */

process deduplicate_valid {
    memory (params.large ? '400 GB' : '64 GB')
    queue = params.longqueue

    input:
        path("valid.genbank.json")
        path("verified.genbank.json")

    output:
        path("deduplicated.valid.json"), emit: deduplicated
        path("all.valid.json"), emit: all

    script:
        """
        deduplicator.py ./ valid all ${task.cpus}
        mv deduplicated.valid.all.json deduplicated.valid.json
        """
}


process sorfdb {
    publishDir "${params.db}/sorfdb", mode: 'copy', pattern: "*.gz"
    memory (params.large ? '172 GB' : '64 GB')
    cpus (params.threads >= 4 ? 4 : params.threads)
    // queue = params.longqueue

    input:
        path(genbank)
        path(blastx)
        path(swissprot)
        path(uniprot)
        path(smprot)
        path(taxonomy)

    output:
        tuple path("sorfdb.*.tar.gz"), path("sorfdb.*.tsv.gz"), path("sorfdb.*.fna.gz"), path("sorfdb.*.faa.gz")

   script:
        """
        export_to_db.py --genbank $genbank --blastx $blastx --swissprot $swissprot --uniprot $uniprot --smprot $smprot  \
        --taxonomy $taxonomy --pfam ${params.db}/pfam/pfam_short.h3m --threads ${task.cpus}
        """
}


process sorfdb_taxonomy {
    publishDir "${params.db}/sorfdb", mode: 'copy', pattern: "taxonomy.*.tsv"
    publishDir "${params.db}/sorfdb", mode: 'copy', pattern: "*.html"
    conda "$baseDir/envs/db-taxonomy.yml"

    input:
        tuple path(tar), path(tsv), path(fna), path(faa)

    output:
        path("taxonomy.*.tsv"), emit: taxtsv
        path("*.html"), emit: krona

   script:
        """
        pigz -dck $tsv | tail -n +2 | cut -f 8-12,18 | sort > tmp_taxonomy.tsv
        taxonomy2krona.py --input tmp_taxonomy.tsv
        ktImportText taxonomy.total.tsv -o sorfdb.taxonomy.total.html
        ktImportText taxonomy.non-redundant.tsv -o sorfdb.taxonomy.non-redundant.total.html
        """
}


process blast_export {
    publishDir "${params.db}/sorfdb", mode: 'copy', pattern: "*.gz"
    memory { 32.GB * task.attempt }
    maxRetries 4

    input:
        path(genbank)
        path(blastx)
        path(swissprot)
        path(uniprot)
        path(smprot)

    output:
        path('sprot.clustering.faa.gz')

   script:
        """
        export_for_blast.py --genbank $genbank --blastx $blastx --swissprot $swissprot --uniprot $uniprot \
        --smprot $smprot --threads ${task.cpus}
        """
}

// EOF
