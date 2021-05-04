nextflow.enable.dsl = 2

include { FASTQC as FASTQC_UNTRIMMED } from "./modules/fastqc/fastqc.nf" addParams(resultsDir: "${params.outdir}/fastqc_untrimmed")
include { FASTQC as FASTQC_TRIMMED } from "./modules/fastqc/fastqc.nf" addParams(resultsDir: "${params.outdir}/fastqc_trimmed")
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
    gatk38_jar_ch = Channel.value(java.nio.file.Paths.get("$params.gatk38_jar"))
    env_user_ch = Channel.value("root")

   FASTQC_UNTRIMMED(reads_ch) // DONE
    TRIMMOMATIC(reads_ch) // DONE
//    FASTQC_TRIMMED(TRIMMOMATIC.out)

    MTBSEQ_PER_SAMPLE(TRIMMOMATIC.out,
            gatk38_jar_ch,
            env_user_ch) // TODO


    samples_tsv_file_ch = MTBSEQ_PER_SAMPLE.out[0]
            .collect()
            .flatten().map { n ->  "$n" + "\t" + "$params.mtbseq_library_name" + "\n"  }
            .collectFile(name: 'samples.tsv', newLine: false, storeDir: "$params.resultsDir_mtbseq_cohort")

    MTBSEQ_COHORT(
            samples_tsv_file_ch,
            MTBSEQ_PER_SAMPLE.out[2].collect(),
            MTBSEQ_PER_SAMPLE.out[3].collect(),
            gatk38_jar_ch,
            env_user_ch,
    ) // TODO


    RDANALYZER(TRIMMOMATIC.out)
    SPOTYPING(TRIMMOMATIC.out)
    SPADES(TRIMMOMATIC.out)
    PROKKA(SPADES.out.prokka_input)
//    QUAST(SPADES.out.quast_input.collect()) // TODO

    TBPROFILER_PROFILE(TRIMMOMATIC.out)

   TBPROFILER_COLLATE(TBPROFILER_PROFILE.out.collect())


}


workflow SPADES_PROKKA_WF {

    reads_ch = Channel.fromFilePairs(params.reads)

    TRIMMOMATIC(reads_ch)
    SPADES(TRIMMOMATIC.out)
    PROKKA(SPADES.out)

}




workflow test {
    reads_ch = Channel.fromFilePairs(params.reads)

    reads_ch.view()
}

