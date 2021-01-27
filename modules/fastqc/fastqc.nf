nextflow.enable.dsl = 2

params.resultsDir = "${params.outdir}/fastqc"
params.saveMode = 'copy'
params.shouldPublish = true

process FASTQC {
    tag "${genomeName}"
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish
    container 'quay.io/biocontainers/fastqc:0.11.9--0'

    input:
    tuple val(genomeName), path(genomeReads)

    output:
    tuple path('*.html'), path('*.zip')


    script:

    """
    fastqc *fastq*
    """

    stub:
    """
    touch ${genomeName}.html
    touch ${genomeName}.zip
    """
}




