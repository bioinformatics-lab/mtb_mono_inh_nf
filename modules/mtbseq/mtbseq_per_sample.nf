nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar

params.results_dir = "${params.outdir}/mtbseq/per_sample"
params.save_mode = 'copy'
params.should_publish = true

process MTBSEQ_PER_SAMPLE {
    tag "${genomeFileName}"
    publishDir params.results_dir, pattern: "${genomeFileName}_results", mode: params.save_mode, enabled: params.should_publish
    // TODO port to errorStrategy and maxRetries
    validExitStatus 0, 1, 2


    input:
    tuple val(genomeFileName), path("${genomeFileName}_${params.mtbseq_library_name}_R?.fastq.gz")
    path(gatk_jar)
    env USER

    output:
    val("${genomeFileName}") // Genome name
    path("""${genomeFileName}_results""") // Folder
    path("""${genomeFileName}_results/Called/*tab""")
    path("""${genomeFileName}_results/Position_Tables/*tab""")

    script:
    """
    set +e

    gatk-register ${gatk_jar}

    mkdir ${genomeFileName}_results

    MTBseq --step TBfull --thread ${task.cpus}

    mv  Amend ./${genomeFileName}_results/
    mv  Bam ./${genomeFileName}_results/
    mv  Called ./${genomeFileName}_results/
    mv  Classification ./${genomeFileName}_results/
    mv  GATK_Bam ./${genomeFileName}_results/
    mv  Groups ./${genomeFileName}_results/
    mv  Joint ./${genomeFileName}_results/
    mv  Mpileup ./${genomeFileName}_results/
    mv  Position_Tables ./${genomeFileName}_results/
    mv  Statistics ./${genomeFileName}_results/
    """

    stub:
    """
    echo "MTBseq --step TBfull --thread ${task.cpus}"
    mkdir ${genomeFileName}_results
    mkdir ${genomeFileName}_results/Called -p
    touch ${genomeFileName}_results/Called/${genomeFileName}_somelib.gatk_position_uncovered_cf4_cr4_fr75_ph4_outmode000.tab
    touch ${genomeFileName}_results/Called/${genomeFileName}_somelib.gatk_position_variants_cf4_cr4_fr75_ph4_outmode000.tab

    mkdir ${genomeFileName}_results/Position_Tables -p
    touch ${genomeFileName}_results/Position_Tables/${genomeFileName}_somelib.gatk_position_table.tab
    """

}

workflow test {


    include { TRIMMOMATIC } from "../../modules/trimmomatic/trimmomatic.nf"
    include { MTBSEQ_COHORT } from "../../modules/mtbseq/mtbseq_cohort.nf"

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
