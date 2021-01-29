nextflow.enable.dsl=2


params.saveMode = 'copy'
params.resultsDir = "${params.outdir}/quast"
params.shouldPublish = true


process QUAST {
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish
    container 'quay.io/biocontainers/quast:5.0.2--py37pl526hb5aa323_2'


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

    mkdir basic_stats/
    mkdir icarus_viewers/
    
    touch icarus.html
    touch report.html
    touch report.pdf
    touch report.tex
    touch report.tsv
    touch report.txt
    touch transposed_report.tex
    touch transposed_report.tsv
    touch transposed_report.txt 
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

