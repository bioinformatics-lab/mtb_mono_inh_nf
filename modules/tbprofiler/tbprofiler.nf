params.resultsDir = 'results/tbProfiler'
params.saveMode = 'copy'

params.saveMode = 'copy'
params.resultsDir = "${params.outdir}/tbprofiler"
params.shouldPublish = true


process TBPROFILER_PROFILE {
    publishDir "$params.resultsDir/profile", mode: params.saveMode, enabled: params.shouldPublish
    container 'quay.io/biocontainers/tb-profiler:2.8.6--pypy_0'
    cpus 4
    memory "7 GB"

    input:
    tuple val(genomeName), file(genomeReads)

    output:
    tuple path("""${genomeName}.results.txt"""),
                 path("""${genomeName}.results.json""")


    script:

    """
    tb-profiler profile -1 ${genomeReads[0]} -2 ${genomeReads[1]}  -t 4 -p $genomeName --txt
    cp results/* ./
    """

}


process TBPROFILER_COLLATE {
    publishDir "$params.resultsDir/profile", mode: params.saveMode, enabled: params.shouldPublish
    container 'quay.io/biocontainers/tb-profiler:2.8.6--pypy_0'
    cpus 4
    memory "7 GB"
    input:
    path("""${params.resultsDir}/results""")

    output:
    tuple path("""${params.resultsDir}/tbprofiler.dr.indiv.itol.txt"""),
            path("""${params.resultsDir}/tbprofiler.dr.itol.txt"""),
            path("""${params.resultsDir}/tbprofiler.json"""),
            path("""${params.resultsDir}/tbprofiler.lineage.itol.txt"""),
            path("""${params.resultsDir}/tbprofiler.txt"""),
            path("""${params.resultsDir}/tbprofiler.variants.txt""")


    script:

    """
    cd $params.resultsDir
    tb-profiler update_tbdb
    tb-profiler collate
    """

}

