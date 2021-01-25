nextflow.enable.dsl = 2

//include { FASTQC } from "./modules/fastqc/fastqc.nf"
include { MTBSEQ } from "./modules/mtbseq/mtbseq.nf"
//include { MULTIQC } from "./modules/multiqc/multiqc.nf"
//include { PROKKA } from "./modules/prokka/prokka.nf"
//include { QUAST } from "./modules/quast/quast.nf"
//include { RDANALYZER } from "./modules/rdanalyzer/rdanalyzer.nf"
//include { SPADES } from "./modules/spades/spades.nf"
//include { SPOTYPING } from "./modules/spotyping/spotyping.nf"
//include { TBPROFILER } from "./modules/tbprofiler/tbprofiler.nf"
include { TRIMMOMATIC } from "./modules/trimmomatic/trimmomatic.nf"


workflow mtbseq {
    reads_ch = Channel.fromFilePairs(params.reads)
    gatk38_jar_ch = Channel.value(params.gatk38_jar)
    env_user_ch = Channel.value("root")

    TRIMMOMATIC(reads_ch)
    MTBSEQ(TRIMMOMATIC.out,
            gatk38_jar_ch,
            env_user_ch)

}

