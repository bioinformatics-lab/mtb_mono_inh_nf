nextflow.enable.dsl = 2

include { PROKKA } from "../../modules/prokka/prokka.nf"
include { SPADES } from "../../modules/spades/spades.nf"
include { TRIMMOMATIC } from "../../modules/trimmomatic/trimmomatic.nf"


workflow SPADES_PROKKA_WF {
    take:
    reads_ch = Channel.fromFilePairs(params.reads)

    main:
    TRIMMOMATIC(reads_ch)
    SPADES(TRIMMOMATIC.out)
    PROKKA(SPADES.out.prokka_input)

}
