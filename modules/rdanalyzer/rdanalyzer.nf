
nextflow.enable.dsl= 2

params.saveMode = 'copy'
params.resultsDir = "${params.outdir}/rdanalyzer"
params.shouldPublish = true


process RDANALYZER {
    tag "${genomeFileName}"
    container 'quay.io/bioinformatics_playground/rd_analyzer:0.0.1'
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish


    input:
    tuple val(genomeFileName), path(genomeReads)

    output:
    tuple path("*result"), path("*depth")


    script:

    """
    python  /RD-Analyzer/RD-Analyzer.py  -o ${genomeFileName} ${genomeReads[0]} ${genomeReads[1]}
    """

    stub:
    """
    echo "python /RD-Analyzer/RD-Analyzer.py -o ${genomeFileName} ${genomeReads[0]} ${genomeReads[1]}"

    touch ${genomeFileName}.result
    touch ${genomeFileName}.depth
    """

}


workflow test {

include { TRIMMOMATIC } from "../trimmomatic/trimmomatic.nf"

input_ch = Channel.fromFilePairs("$launchDir/test_data/*_{1,2}.fastq.gz")

TRIMMOMATIC(input_ch)

RD_ANALYZER(TRIMMOMATIC.out)



}
