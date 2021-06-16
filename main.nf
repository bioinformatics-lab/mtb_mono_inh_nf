nextflow.enable.dsl = 2

include { FASTQC as FASTQC_UNTRIMMED } from "./modules/fastqc/fastqc.nf" addParams(resultsDir: "${params.outdir}/fastqc_untrimmed")
include { FASTQC as FASTQC_TRIMMED } from "./modules/fastqc/fastqc.nf" addParams(resultsDir: "${params.outdir}/fastqc_trimmed")
include { MULTIQC as MULTIQC_TRIMMED } from "./modules/multiqc/multiqc.nf" addParams(resultsDir: "${params.outdir}/multiqc_trimmed", fastqcResultsDir: "${params.outdir}/fastqc_trimmed")
include { MULTIQC as MULTIQC_UNTRIMMED } from "./modules/multiqc/multiqc.nf" addParams(resultsDir: "${params.outdir}/multiqc_untrimmed", fastqcResultsDir: "${params.outdir}/fastqc_untrimmed")
include { MTBSEQ_PER_SAMPLE; MTBSEQ_COHORT } from "./modules/mtbseq/mtbseq.nf"
include { MULTIQC } from "./modules/multiqc/multiqc.nf"
include { PROKKA } from "./modules/prokka/prokka.nf"
include { RDANALYZER } from "./modules/rdanalyzer/rdanalyzer.nf"
include { SPADES } from "./modules/spades/spades.nf"
include { SPOTYPING } from "./modules/spotyping/spotyping.nf"
include { QUAST } from "./modules/quast/quast.nf"
include { TBPROFILER_PROFILE; TBPROFILER_COLLATE } from "./modules/tbprofiler/tbprofiler.nf"
include { TRIMMOMATIC } from "./modules/trimmomatic/trimmomatic.nf"


workflow {
    reads_ch = Channel.fromFilePairs(params.reads)

    env_user_ch = Channel.value("root")


    FASTQC_UNTRIMMED(reads_ch)
    MULTIQC_UNTRIMMED(FASTQC_UNTRIMMED.out.flatten().collect())
    TRIMMOMATIC(reads_ch)
    FASTQC_TRIMMED(TRIMMOMATIC.out)
    MULTIQC_TRIMMED(FASTQC_TRIMMED.out.flatten().collect())


    MTBSEQ_PER_SAMPLE(TRIMMOMATIC.out,
                      params.gatk38_jar,
                      env_user_ch)


    samples_tsv_file_ch = MTBSEQ_PER_SAMPLE.out[0]
            .collect()
            .flatten().map { n -> "$n" + "\t" + "$params.mtbseq_library_name" + "\n" }
            .collectFile(name: 'samples.tsv', newLine: false, storeDir: "$params.resultsDir_mtbseq_cohort")

    MTBSEQ_COHORT(
            samples_tsv_file_ch,
            MTBSEQ_PER_SAMPLE.out[2].collect(),
            MTBSEQ_PER_SAMPLE.out[3].collect(),
            params.gatk38_jar,
            env_user_ch,
    )

    RDANALYZER(TRIMMOMATIC.out)
    SPOTYPING(TRIMMOMATIC.out)
    SPADES(TRIMMOMATIC.out)
    PROKKA(SPADES.out.prokka_input)

    QUAST(SPADES.out.quast_input.collect())

    TBPROFILER_PROFILE(TRIMMOMATIC.out)
    TBPROFILER_COLLATE(TBPROFILER_PROFILE.out.collect())


}


workflow SAME_PERSON_GENOMES_WF {
    reads_ch = Channel.fromFilePairs(params.reads)

    FASTQC_UNTRIMMED(reads_ch)
    MULTIQC_UNTRIMMED(FASTQC_UNTRIMMED.out.flatten().collect())
    TRIMMOMATIC(reads_ch)
    FASTQC_TRIMMED(TRIMMOMATIC.out)
    MULTIQC_TRIMMED(FASTQC_TRIMMED.out.flatten().collect())

    SPOTYPING(TRIMMOMATIC.out)

    TBPROFILER_PROFILE(TRIMMOMATIC.out)
    TBPROFILER_COLLATE(TBPROFILER_PROFILE.out.collect())


}


workflow AWS_WF {
    reads_ch = Channel.fromFilePairs(params.reads)
    env_user_ch = Channel.value("root")


    // FASTQC_UNTRIMMED(reads_ch)
    // TRIMMOMATIC(reads_ch)
    // FASTQC_TRIMMED(TRIMMOMATIC.out)
    //

    trimmomatic_ch = Channel.fromFilePairs("s3://bucket-eu-west-1-ireland/mono_inh_cristina/results/trimmomatic/*_{R1,R2}.p.fastq.gz")

    MTBSEQ_PER_SAMPLE(trimmomatic_ch,
                      params.gatk38_jar,
                      env_user_ch)

    samples_tsv_file_ch = MTBSEQ_PER_SAMPLE.out[0]
            .collect()
            .flatten().map { n -> "$n" + "\t" + "$params.mtbseq_library_name" + "\n" }
            .collectFile(name: 'samples.tsv', newLine: false, storeDir: "$params.resultsDir_mtbseq_cohort")

    MTBSEQ_COHORT(
            samples_tsv_file_ch,
            MTBSEQ_PER_SAMPLE.out[2].collect(),
            MTBSEQ_PER_SAMPLE.out[3].collect(),
            params.gatk38_jar,
            env_user_ch
    )

    // RDANALYZER(TRIMMOMATIC.out)
    SPOTYPING(trimmomatic_ch)
    SPADES(trimmomatic_ch)
    PROKKA(SPADES.out.prokka_input)

    QUAST(SPADES.out.quast_input.collect())

    TBPROFILER_PROFILE(trimmomatic_ch)
    TBPROFILER_COLLATE(TBPROFILER_PROFILE.out.collect())


}


workflow SPADES_PROKKA_WF {

    reads_ch = Channel.fromFilePairs(params.reads)

    TRIMMOMATIC(reads_ch)
    SPADES(TRIMMOMATIC.out)
    PROKKA(SPADES.out.prokka_input)

}


workflow test_wf {
    trimmomatic_ch = Channel.fromFilePairs("s3://bucket-eu-west-1-ireland/mono_inh_cristina/results/trimmomatic/*_{R1,R2}.p.fastq.gz")

    trimmomatic_ch.view()
}
