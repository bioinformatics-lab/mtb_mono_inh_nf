nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/spotyping"
params.save_mode = 'copy'
params.R2 = false
params.should_publish = true


process SPOTYPING {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path(genomeReads)

    output:
    tuple file('*.txt'), file('SITVIT*.xls')

    script:
    genomeReadToBeAnalyzed = params.R2 ? genomeReads[1] : genomeReads[0]

    """
    python /SpoTyping-v2.0/SpoTyping-v2.0-commandLine/SpoTyping.py ./${genomeReadToBeAnalyzed} -o ${genomeFileName}.txt
    """

    stub:
    genomeReadToBeAnalyzed = params.R2 ? genomeReads[1] : genomeReads[0]
    """
    echo "python /SpoTyping-v2.0/SpoTyping-v2.0-commandLine/SpoTyping.py ./${genomeReadToBeAnalyzed} -o ${genomeFileName}.txt"
    touch ${genomeFileName}.txt
    touch SITVIT_${genomeFileName}.xls
    """

}


workflow test {

    params.TRIMMOMATIC = [
            shouldPublish: false
    ]


    include { TRIMMOMATIC } from "../trimmomatic/trimmomatic.nf" addParams(params.TRIMMOMATIC)

    input_ch = Channel.fromFilePairs("$launchDir/test_data/*_{1,2}.fastq.gz")

    TRIMMOMATIC(input_ch)

    SPOTYPING(TRIMMOMATIC.out)


}
