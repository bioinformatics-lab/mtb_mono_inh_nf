nextflow.enable.dsl = 2

include { FASTQC as FASTQC_UNTRIMMED } from "./modules/fastqc/fastqc.nf" addParams(resultsDir: "${params.outdir}/fastqc_untrimmed")
include { FASTQC as FASTQC_TRIMMED } from "./modules/fastqc/fastqc.nf" addParams(resultsDir: "${params.outdir}/fastqc_trimmed")
include { MTBSEQ } from "./modules/mtbseq/mtbseq.nf"
//include { MULTIQC } from "./modules/multiqc/multiqc.nf"
include { PROKKA } from "./modules/prokka/prokka.nf"
//include { QUAST } from "./modules/quast/quast.nf"
include { RDANALYZER } from "./modules/rdanalyzer/rdanalyzer.nf"
include { SPADES } from "./modules/spades/spades.nf"
include { SPOTYPING } from "./modules/spotyping/spotyping.nf"
//include { TBPROFILER_PROFILE; TBPROFILER_COLLATE } from "./modules/tbprofiler/tbprofiler.nf"
include { TRIMMOMATIC } from "./modules/trimmomatic/trimmomatic.nf"




workflow {
    reads_ch = Channel.fromFilePairs(params.reads)
    gatk38_jar_ch = Channel.value(params.gatk38_jar)
    env_user_ch = Channel.value("root")


    FASTQC_UNTRIMMED(reads_ch)

    TRIMMOMATIC(reads_ch)
    FASTQC_TRIMMED(TRIMMOMATIC.out)
//    MTBSEQ(TRIMMOMATIC.out,
//            gatk38_jar_ch,
//            env_user_ch)

    RDANALYZER(TRIMMOMATIC.out)
    SPOTYPING(TRIMMOMATIC.out)
    SPADES(TRIMMOMATIC.out)
    PROKKA(SPADES.out)

//    TBPROFILER_PROFILE(TRIMMOMATIC.out)
//    TBPROFILER_COLLATE(TBPROFILER_PROFILE.out)


}



workflow mtbseq {
    reads_ch = Channel.fromFilePairs(params.reads)
    gatk38_jar_ch = Channel.value(params.gatk38_jar)
    env_user_ch = Channel.value("root")

    TRIMMOMATIC(reads_ch)
    MTBSEQ(TRIMMOMATIC.out,
            gatk38_jar_ch,
            env_user_ch)

}



workflow test {
    reads_ch = Channel.fromFilePairs(params.reads)

    reads_ch.view()
}

