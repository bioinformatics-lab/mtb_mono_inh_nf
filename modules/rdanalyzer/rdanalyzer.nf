
nextflow.enable.dsl= 2

params.saveMode = 'copy'
params.resultsDir = "${params.outdir}/rd_analyzer"
params.shouldPublish = true


process RDANALYZER {
    tag "${genomeFileName}"
    container 'quay.io/bioinformatics_playground/rd_analyzer:0.0.1'
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish
    cpus 8
    memory "16 GB"


    input:
    tuple val(genomeFileName), path(genomeReads)

    output:
    tuple path("*result"), path("*depth")


    script:
    genomeName = genomeFileName.toString().split("\\_")[0]

    """
    python  /RD-Analyzer/RD-Analyzer.py  -o ${genomeName} ${genomeReads[0]} ${genomeReads[1]}
    """
}


workflow test {

include { TRIMMOMATIC } from "../trimmomatic/trimmomatic.nf"

input_ch = Channel.fromFilePairs("$launchDir/test_data/*_{1,2}.fastq.gz")

TRIMMOMATIC(input_ch)

RD_ANALYZER(TRIMMOMATIC.out)



}
