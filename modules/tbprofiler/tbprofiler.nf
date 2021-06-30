nextflow.enable.dsl = 2

params.results_dir_tbprofiler_per_sample = "${params.outdir}/tbprofiler/samples"
params.save_mode_tbprofiler_per_sample = 'copy'
params.should_publish_tbprofiler_per_sample = true


process TBPROFILER_PROFILE {
    tag "${genomeName}"
    publishDir params.results_dir_tbprofiler_per_sample, mode: params.save_mode_tbprofiler_per_sample, enabled: params.should_publish_tbprofiler_per_sample

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
    echo "tb-profiler profile -1 ${genomeReads[0]} -2 ${genomeReads[1]}  -t ${task.cpus} -p $genomeName --txt"

    mkdir results
    touch results/"${genomeName}.results.txt"
    touch results/"${genomeName}.results.json"
    """

}


params.results_dir_tbprofiler_cohort = "${params.outdir}/tbprofiler/cohort"
params.save_mode_tbprofiler_cohort = 'copy'
params.should_publish_tbprofiler_cohort = true


process TBPROFILER_COLLATE {
    publishDir params.results_dir_tbprofiler_cohort, mode: params.save_mode_tbprofiler_cohort, enabled: params.should_publish_tbprofiler_cohort

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
