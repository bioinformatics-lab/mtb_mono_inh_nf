nextflow.enable.dsl= 2

params.resultsDir_tbprofiler_per_sample = "${params.outdir}/tbprofiler/samples"
params.saveMode_tbprofiler_per_sample = 'copy'
params.shouldPublish_tbprofiler_per_sample = true


process TBPROFILER_PROFILE {
    publishDir params.resultsDir_tbprofiler_per_sample, mode: params.saveMode_tbprofiler_per_sample, enabled: params.shouldPublish_tbprofiler_per_sample
    container 'quay.io/biocontainers/tb-profiler:2.8.6--pypy_0'

    input:
    tuple val(genomeName), path(genomeReads)

    output:
    tuple path("results/*txt"), path("results/*json")


    script:
    """
    tb-profiler profile -1 ${genomeReads[0]} -2 ${genomeReads[1]}  -t ${task.cpus} -p $genomeName --txt
    """

    stub:
    """
    mkdir results
    touch results/"${genomeName}.results.txt"
    touch results/"${genomeName}.results.json"
    """

}


params.resultsDir_tbprofiler_cohort = "${params.outdir}/tbprofiler/cohort"
params.saveMode_tbprofiler_cohort = 'copy'
params.shouldPublish_tbprofiler_cohort = true


process TBPROFILER_COLLATE {
    publishDir params.resultsDir_tbprofiler_cohort, mode: params.saveMode_tbprofiler_cohort, enabled: params.shouldPublish_tbprofiler_cohort
    container 'quay.io/biocontainers/tb-profiler:2.8.6--pypy_0'

    input:
    path("results/*")

    output:
    path("tbprofiler*")

    script:
    """
    tb-profiler update_tbdb
    tb-profiler collate
    cp tbprofiler.txt tbprofiler_cohort_report.tsv
    """

    stub:
    """
    touch tbprofiler_cohort_report.tsv
    """
}


workflow test {
    include { TRIMMOMATIC } from "../../modules/trimmomatic/trimmomatic.nf"

    input_reads_ch = Channel.fromFilePairs("$launchDir/data/mock_data/*_{R1,R2}*fastq.gz")

    TRIMMOMATIC(input_reads_ch)

    TBPROFILER_PROFILE(TRIMMOMATIC.out)

    TBPROFILER_COLLATE(
            TBPROFILER_PROFILE.out.collect()
    )


}
