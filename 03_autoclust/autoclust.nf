#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import java.nio.file.*

params.selector = 'eff'
params.srv = 30
params.id = params.srv
params.cov = params.srv
params.qcov = params.cov
params.scov = params.cov
params.minsize = 5
params.block = 100.KB

println("Executor=" + params.executor)
println("Threads=" + params.threads)
println("SRV=" + params.srv)
params.mamba ? println("Using Mamba.") : println("Using Conda.")
println("\n")


workflow {
    fasta = Channel.fromPath(params.input)
    fasta_chunks = Channel.fromPath(params.input).splitFasta(size: params.block, file: true)
    sorfdb_tsv = Channel.fromPath("${params.db}/sorfdb/sorfdb.*.tsv.gz")

    inflation = Channel.of(
            1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
            3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0
          )

    count_fasta_entries(params.input)
    fasta_count = count_fasta_entries.out.map{f -> f.text.toInteger()}

    makeblastdb(fasta)
    blast_input = fasta_chunks.combine(makeblastdb.out)
    blastp(blast_input, fasta_count)
    srv_transform(blastp.out)

    // save raw blast output
    compress(blastp.out.collectFile(name: 'blast.tsv',
                                    keepHeader: false,
                                    cache: 'lenient',
                                    sort: 'hash'))

    // prune singletons, pairs and nodes that are connected only to one node of a graph
    prune(srv_transform.out.collectFile(name: 'srv.tsv',
                                        keepHeader: false,
                                        cache: 'lenient',
                                        sort: 'hash'))  // storeDir: "${params.db}/sorfdb",

    // load to binary mcl matrix; make the matrix symmetric (mean of SRVs)
    mcxload(prune.out)
    // select k best edges of a node without introducing singletons
    mcxalter(mcxload.out)
    // write the pruned, symmetric and knn reduced data back to a abc TSV file
    mcxdump_abc(mcxalter.out)
    // separate unlinked subgraphs with similar attributes into separate files
    separate(mcxdump_abc.out)

    // load separated abc TSVs to binary mcl matrices
    // input: id_subgraph, abc TSV
    // ouput: id_subgraph, mcx, tab
    mcxload2(separate.out.flatten().map{[it.getBaseName().tokenize('.')[0].toInteger(), it]})

    // mcl clustering for all separated graphs and inflations
    // input: inflation, id_subgraph, mcx, tab
    // output: id_subgraph, inflation, mci.I
    mcl(inflation.combine(mcxload2.out))

    // look up the optimal inflation for the separated subgrpahs with clm info based on modularity
    // input: id_subgraph, []inflation, ...], [mci.I, ...], mcx, tab
    // output: id_subgraph, (best) inflation, (best) id_subgraph.mci.I, tab
    clm(mcl.out.groupTuple().join(mcxload2.out))

    // dump the clusters for the separated subgrpahs
    // input: id_subgraph, inflation, mci.I, tab
    // output: id_subgraph, inflation, cluster
    mcxdump_clusters(clm.out.cluster)

    // collect all clusters of the separated subgrpahs into one file
    // input: id_subgraph, inflation, clusters -> clusters
    // output: all clusters in one file
    autoclust = mcxdump_clusters.out.map{it[2]}.collectFile(name: 'auto.cluster.tsv',
                                                            storeDir: "${params.db}/sorfdb/raw/",
                                                            keepHeader: false,
                                                            cache: 'lenient',
                                                            sort: 'hash')

    // group mci.I and tab files together by their subgraph id for all inflation values
    // input: id_subgraph, inflation, mci.I; id_subgraph, mcx, tab
    // output: inflation, mci.I, tab
    clusters_by_inflation = mcl.out.groupTuple().join(mcxload2.out).map{
        items -> {
        def combined = []
            if (items[1] instanceof List) {
                items[1].eachWithIndex{it, i -> {  //
                    combined.add([it,
                                  items[2][i],
                                  items[4]])}
                }
            } else {
                combined.add([items[1],
                              items[2],
                              items[4]])
            }
            return combined
        }
    }.flatten().collate(3)

    // dump clusters for all inflation values for each subgraph
    // input: inflation, mci.I, tab
    // output: inflation, clusters
    mcxdump_for_cytoscape(clusters_by_inflation)

    // group clusters by inflation value
    clusters_by_inflation = mcxdump_for_cytoscape.out.groupTuple().map{it[1]}.flatten().collectFile(
        keepHeader: false,
        cache: 'lenient',
        sort: 'hash'
    ).collect()

    // export the graph with cluster information for cytoscape
    // input: srv.tsv; sorfdb; autoclust_clusters; clusters_for_all_inflations
    // output: grraph_with_cluster_information
    cytoscape(mcxdump_abc.out,
              sorfdb_tsv,
              autoclust,
              clusters_by_inflation)

    // export all clusters to fasta files
    // input: id_subgraph, inflation, clusters
    // output: id_subgraph.cluster_size.counter.fasta
    export_clusters_to_fasta(mcxdump_clusters.out)


    // add subgraph id and the number of sequences for each fasta file
    // input: id_subgraph.cluster_size.counter.fasta
    // output: id_subgraph.cluster_size.counter, sequence_count, id_subgraph.cluster_size.counter.fasta
    cluster_fasta = export_clusters_to_fasta.out.map{items -> {
                                                        def combined = []
                                                        if (items instanceof List) {
                                                            items.each{combined.add([it.getBaseName(),  // id_subgraph.cluster_size.counter
                                                                                     it.getBaseName().tokenize('.')[1].toInteger(),  // size
                                                                                     it])}  // path
                                                        } else {
                                                            combined.add([items.getBaseName(),  // id_subgraph.cluster_size.counter
                                                                          items.getBaseName().tokenize('.')[1].toInteger(),  // size
                                                                          items])  // path
                                                        }
                                                        return combined
                                                    }}.flatten().collate(3)

    // sort the fasta to the according muscle algorithm depending on their number of sequences
    aln = cluster_fasta.branch{
        ppp: it[1] <= 1000
        super5: it[1] > 1000
    }
    
    muscle_ppp(aln.ppp)
    muscle_super5(aln.super5)
    
    
    // build HMMs based on the multiple alignments
    // input: id_subgraph.cluster_size.counter, sequence_count, id_subgraph.cluster_size.counter.fasta, sprot.fasta
    // output: aln.id_subgraph.cluster_size.counter.afa
    hmm_build(muscle_ppp.out.concat(muscle_super5.out).collect(), sorfdb_tsv, fasta)
    compress_hmm(hmm_build.out.hmm)
}


process count_fasta_entries {
    conda "$baseDir/envs/db-clustering.yml"

    input:
        path(faa)

    output:
        path('faa.count')

    script:
    """
    zgrep -c '>' $faa > faa.count
    """
}


process makeblastdb {
    conda "$baseDir/envs/db-clustering.yml"
    memory { 4.GB * task.attempt }

    input:
        path(faa)

    output:
        path('blastdb/')

    script:
    """
    gzip -dck $faa > sprot.faa
    makeblastdb -dbtype prot -hash_index -in sprot.faa -out blastdb
    mkdir blastdb
    mv blastdb.* blastdb/
    """
}


process blastp {
    conda "$baseDir/envs/db-clustering.yml"
    memory { 8.GB * task.attempt }
    maxRetries 3

    input:
        tuple path(faa), path(db)
        val(count)

    output:
        path('blast.tsv')

    script:
    """
    blastp -task blastp-short \
    -query $faa \
    -db $db/blastdb \
    -out blast.tsv \
    -max_target_seqs $count \
    -matrix BLOSUM62 \
    -comp_based_stats 0 \
    -soft_masking false \
    -seg no \
    -evalue 1 \
    -outfmt '6 qseqid sseqid pident length bitscore evalue' \
    -num_threads ${task.cpus}
    """
}


process compress {
    publishDir "${params.db}/sorfdb/raw", mode: 'copy', pattern: "*.gz"
    conda "$baseDir/envs/db-clustering.yml"
    cpus (params.threads >= 10 ? 10: params.threads)
    memory { 1.GB * task.attempt }

    input:
        path('blast.tsv')

    output:
        path("*.gz")

    script:
    """
    pigz -f -p ${task.cpus} -9 blast.tsv
    """
}


process srv_transform {
    conda "$baseDir/envs/db-clustering.yml"
    cpus (params.threads >= 2 ? 2 : params.threads)
    memory { 10.GB * task.attempt }
    maxRetries 3

    input:
        path('blast.tsv')

    output:
        path('srv.tsv')

    script:
    """
    blast_to_srv_matrix.py --input blast.tsv --mode blast --minsrv ${params.srv} --id ${params.id} \
    --query_cover ${params.qcov} --subject_cover ${params.scov} --threads ${task.cpus}
    """
}

process prune {
    conda "$baseDir/envs/db-clustering.yml"

    input:
        path('srv.tsv')

    output:
        path('pruned.srv.tsv')

    script:
    """
    srv_pruning.py --input srv.tsv --low-memory --threads ${task.cpus}
    """
}


process mcxload {
    conda "$baseDir/envs/db-clustering.yml"
    memory { 32.GB * task.attempt }

    input:
        path('srv.tsv')

    output:
        tuple path('srv.mcx'), path('srv.tab')

    script:
    """
    mcxload -ri add -tf 'scale(2),gt(${params.srv})' --write-binary -abc srv.tsv -write-tab srv.tab -o srv.mcx
    """
}


process mcxalter {
    conda "$baseDir/envs/db-clustering.yml"
    memory { 32.GB * task.attempt }

    input:
        tuple path('srv.unaltered.mcx'), path('srv.tab')

    output:
        tuple path('srv.mcx'), path('srv.tab')

    script:
    """
    TOP=\$(mcx query -imx srv.unaltered.mcx -vary-knn 0/1000/100 -t ${task.cpus} | tee top.txt | \
    awk '{if (\$4 ~ /^[0-9]+\$/ && \$4>=1) {print x; exit 0};x=\$0}' | head -n1 | awk '{print \$NF}')
    if [[ \$TOP < 100 ]]; then
        TOP=100
        ln -s srv.unaltered.mcx srv.mcx
    else
        echo "\$TOP" > top.selected.txt
        K=\$(mcx query -imx srv.unaltered.mcx -vary-knn \$((\$TOP-100))/\$TOP/10 -t ${task.cpus} | \
        tee k.txt | awk '{if (\$4 ~ /^[0-9]+\$/ && \$4>=1) {print x; exit 0};x=\$0}' | head -n1 | awk '{print \$NF}')
        echo "\$K" > k.selected.txt

        mcx alter -imx srv.unaltered.mcx -tf "#knn(\$K)" --write-binary -o srv.mcx
    fi
    """
}


process mcxdump_abc {
    publishDir "${params.db}/sorfdb/raw", mode: 'copy', pattern: "pruned.srv.tsv"
    conda "$baseDir/envs/db-clustering.yml"

    input:
        tuple path('srv.mcx'), path('srv.tab')

    output:
        path('pruned.srv.tsv')

    script:
    """
    mcxdump -imx srv.mcx -tab srv.tab -o pruned.srv.tsv
    """
}


process separate {
    conda "$baseDir/envs/db-clustering.yml"
    memory { 16.GB * task.attempt }

    input:
        path('srv.tsv')

    output:
        path('*.subgraph.tsv')

    script:
    """
    separate_graph.py --abc srv.tsv --threads ${task.cpus}
    """
}


process mcxload2 {
    conda "$baseDir/envs/db-clustering.yml"
    memory { 8.GB * task.attempt }

    input:
        tuple val(id), path('srv.tsv')

    output:
        tuple val(id), path('srv.mcx'), path('srv.tab')

    script:
    """
    mcxload --write-binary -abc srv.tsv -write-tab srv.tab -o srv.mcx
    """
}


process mcl {
    conda "$baseDir/envs/db-clustering.yml"
    cpus { params.threads >= 4 * task.attempt * task.attempt ? 4 * task.attempt * task.attempt : params.threads }
    memory { 1.GB * task.attempt }

    input:
        tuple val(inflation), val(id), path('srv.mcx'), path('srv.tab')

    output:
        tuple val(id), val(inflation), path("*.clusters.mci.I*")

    script:
    """
    I=\$(echo $inflation | sed 's/\\.//g')
    echo \$I

    mcl srv.mcx -o ${id}.clusters.mci.I\$I -scheme 7 -I $inflation -te ${task.cpus}
    """
}


process clm {
    publishDir "${params.db}/sorfdb/raw/clm/", mode: 'copy', pattern: "*.txt"
    conda "$baseDir/envs/db-clustering.yml"
    memory { 4.GB * task.attempt }

    input:
        tuple val(id), val(inflations), path(clusters), path('srv.mcx'), path('srv.tab')

    output:
        tuple val(id), env(INFLATION), path("${id}.mci"), path('srv.tab'), emit: cluster
        path("*.txt"), emit: log

    script:
    """
    clm info -o ${id}.clminfo.txt 'srv.mcx' $clusters
    INFLATION=\$(mcl_inflation.py --input ${id}.clminfo.txt --id $id --metric ${params.selector})
    """
}


process mcxdump_clusters {
    conda "$baseDir/envs/db-clustering.yml"

    input:
        tuple val(id), val(inflation), path(mci), path('srv.tab')

    output:
        tuple val(id), val(inflation), path("${id}.cluster.tsv")

    script:
    """
    mcxdump -icl $mci -tabr srv.tab -o ${id}.cluster.tsv
    """
}


process mcxdump_for_cytoscape {
    conda "$baseDir/envs/db-clustering.yml"

    input:
        tuple val(inflation), path(mci), path('srv.tab')

    output:
        tuple val(inflation), path("cluster.${inflation}.tsv")

    script:
    """
    mcxdump -icl $mci -tabr srv.tab -o cluster.${inflation}.tsv
    """
}


process cytoscape {
    publishDir "${params.db}/sorfdb/raw", mode: 'copy', pattern: "srv.clustered.tsv"
    conda "$baseDir/envs/db-clustering.yml"
    memory { 64.GB * task.attempt }

    input:
        path('srv.tsv')
        path(sorfdb)
        path(autoclust)
        path(clusters)

    output:
        path('srv.clustered.tsv')

    script:
    """
    export_to_cytoscape.py --abc srv.tsv --sorfdb $sorfdb --autoclust $autoclust --cluster $clusters
    """
}


process export_clusters_to_fasta {
    conda "$baseDir/envs/db-clustering.yml"

    input:
        tuple val(id), val(inflation), path(clusters)

    output:
        path('*.faa')

    script:
    """
    split_clusters_for_alignment.py --input $clusters --id $id --size ${params.minsize}
    """
}


process muscle_ppp {
    publishDir "${params.db}/sorfdb/raw/aln${params.cov}", mode: 'copy', pattern: "aln.*.afa"
    conda "$baseDir/envs/db-clustering.yml"
    cpus (params.threads >= 8 ? 8 : params.threads)
    memory { 2.GB * task.attempt }
    maxRetries = 5

    input:
        tuple val(identifier), val(size), path(fasta)

    output:
        path("aln.${identifier}.afa")

    script:
    """
    muscle -amino -align $fasta -output aln.${identifier}.afa -threads ${task.cpus}
    """
}


process muscle_super5 {
    publishDir "${params.db}/sorfdb/raw/aln${params.cov}", mode: 'copy', pattern: "aln.*.afa"
    conda "$baseDir/envs/db-clustering.yml"
    cpus (params.threads >= 16 ? 16 : params.threads)
    memory { 4.GB * task.attempt }
    maxRetries = 3

    input:
        tuple val(identifier), val(size), path(fasta)

    output:
        path("aln.${identifier}.afa")

    script:
    """
    muscle -amino -super5 $fasta -output aln.${identifier}.afa -threads ${task.cpus}
    """
}


process hmm_build {
    publishDir "${params.db}/sorfdb", mode: 'copy', pattern: "sorfdb.clusters.*.tsv.gz"
    conda "$baseDir/envs/db-clustering.yml"
    memory { 64.GB * task.attempt }

    input:
        path(aln)
        path(sorfdb)
        path(fasta)

    output:
        path("*.hmm"), emit: hmm
        path("sorfdb.clusters.*.tsv.gz"), emit: clusters

    script:
    """
    alignement_to_hmm.py --sorfdb $sorfdb --proteins $fasta --alignments ./
    """
}

process compress_hmm {
    publishDir "${params.db}/sorfdb", mode: 'copy', pattern: "*.hmm.gz"
    conda "$baseDir/envs/db-clustering.yml"
    memory { 8.GB * task.attempt }
    cpus (params.threads >= 4 ? 4 : params.threads)

    input:
        path(hmm)

    output:
        path("*.hmm.gz")

    script:
    """
    pigz -f -p ${task.cpus} -9 $hmm
    """
}
