nextflow.enable.dsl=2


params.resultsDir = "${params.outdir}/spades"
params.shouldPublish= true
params.saveMode = 'copy'

process SPADES {
    tag "${genomeName}"
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish

    input:
    tuple val(genomeName), path(genomeReads)

    output:
    path "${genomeName}_contigs.fasta", emit: quast_input
    tuple val(genomeName), path("${genomeName}_contigs.fasta"), emit: prokka_input


    script:

    """
    spades.py -k 21,33,55,77 --careful --only-assembler --pe1-1 ${genomeReads[0]} --pe1-2 ${genomeReads[1]} -o ${genomeName} -t ${task.cpus}
    cp ${genomeName}/contigs.fasta ${genomeName}_contigs.fasta 
    """

    stub:
    """
    echo "spades.py -k 21,33,55,77 --careful --only-assembler --pe1-1 ${genomeReads[0]} --pe1-2 ${genomeReads[1]} -o ${genomeName} -t ${task.cpus}"

    touch ${genomeName}_scaffolds.fasta
    touch ${genomeName}_contigs.fasta 
    """
}



workflow test {


include { TRIMMOMATIC } from "$launchDir/modules/trimmomatic/trimmomatic.nf"

input_reads_ch = Channel.fromFilePairs("$launchDir/data/mock_data/*_{R1,R2}*fastq.gz")

TRIMMOMATIC(input_reads_ch)

SPADES(TRIMMOMATIC.out)

SPADES.out.collect().view()
}
