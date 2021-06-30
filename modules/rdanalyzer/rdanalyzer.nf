nextflow.enable.dsl = 2

params.save_mode = 'copy'
params.results_dir = "${params.outdir}/rdanalyzer"
params.should_publish = true


process RDANALYZER {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


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
