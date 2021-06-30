nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar

params.results_dir = "${params.outdir}/mtbseq/cohort"
params.save_mode = 'copy'
params.should_publish = true

process MTBSEQ_COHORT {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish
    // TODO port to errorStrategy and maxRetries
    validExitStatus 0, 1, 2

    input:
    path(samples_tsv_ch)
    path("Called/*")
    path("Position_Tables/*")
    path(gatk_jar)
    env USER

    output:
    tuple path("Joint"), path("Amend"), path("Groups")

    script:

    """
    set +e

    gatk-register ${gatk_jar}

    sleep 10

    mkdir Joint
    MTBseq --step TBjoin --samples ${samples_tsv_ch} --project ${params.mtbseq_project_name}

    mkdir Amend
    MTBseq --step TBamend --samples ${samples_tsv_ch} --project ${params.mtbseq_project_name}

    mkdir Groups
    MTBseq --step TBgroups --samples ${samples_tsv_ch} --project ${params.mtbseq_project_name}

    """



    stub:

    """
    mkdir Joint Amend Groups
    """

}

workflow test {


    include { TRIMMOMATIC } from "../../modules/trimmomatic/trimmomatic.nf"
    include { MTBSEQ_PER_SAMPLE } from "../../modules/mtbseq/mtbseq_per_sample.nf"

    input_reads_ch = Channel.fromFilePairs("$launchDir/data/mock_data/*_{R1,R2}*fastq.gz")
    gatk_jar_ch = Channel.value("$launchDir/resources/GenomeAnalysisTK.jar")
    user_ch = Channel.value("root")

    TRIMMOMATIC(input_reads_ch)
    MTBSEQ_PER_SAMPLE(TRIMMOMATIC.out,
            gatk_jar_ch,
            user_ch)

//    samples_tsv_file_ch = MTBSEQ_PER_SAMPLE.out[1].collect().flatten().map { n ->  "$n" + "\t" + "$params.mtbseq_library_name" + "\n"  }
//    samples_tsv_file_ch.count().view()
//    samples_tsv_file_ch.view()

    samples_tsv_file_ch = MTBSEQ_PER_SAMPLE.out[0]
            .collect()
            .flatten().map { n -> "$n" + "\t" + "$params.mtbseq_library_name" + "\n" }
            .collectFile(name: 'samples.tsv', newLine: false, storeDir: "$params.results_dir_mtbseq_cohort")


    MTBSEQ_COHORT(
            samples_tsv_file_ch,
            MTBSEQ_PER_SAMPLE.out[2].collect(),
            MTBSEQ_PER_SAMPLE.out[3].collect(),
            gatk_jar_ch,
            user_ch
    )

}
