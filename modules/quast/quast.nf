nextflow.enable.dsl = 2


params.saveMode = 'copy'
params.resultsDir = "${params.outdir}/quast"
params.shouldPublish = true


process QUAST {
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish


    input:
    path(scaffoldFiles)

    output:
    path("quast_results")


    script:

    """
    quast ${scaffoldFiles}
    """

    stub:
    """
    echo "quast ${scaffoldFiles}"

    mkdir quast_results/

    """
}

workflow test {
    include { TRIMMOMATIC } from "../trimmomatic/trimmomatic.nf"
    include { SPADES } from "../spades/spades.nf"

    input_reads_ch = Channel.fromFilePairs("$launchDir/data/mock_data/*_{R1,R2}*fastq.gz")

    TRIMMOMATIC(input_reads_ch)

    SPADES(TRIMMOMATIC.out)

    QUAST(SPADES.out.collect())

}

