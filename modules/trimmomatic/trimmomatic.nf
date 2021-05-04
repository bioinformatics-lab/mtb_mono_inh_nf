nextflow.enable.dsl= 2


params.saveMode = 'copy'
params.resultsDir = "${params.outdir}/trimmomatic"
params.shouldPublish = true


process TRIMMOMATIC {
    tag "${genomeName}"
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish
    container 'quay.io/biocontainers/trimmomatic:0.35--6'

    input:
    tuple val(genomeName), path(genomeReads)

    output:
    tuple val(genomeFileName), path("*_{R1,R2}.p.fastq.gz")

    script:

    genomeFileName = genomeName.toString().split("\\_")[0]

    fq_1_paired = genomeFileName + '_R1.p.fastq.gz'
    fq_1_unpaired = genomeFileName + '_R1.s.fastq.gz'
    fq_2_paired = genomeFileName + '_R2.p.fastq.gz'
    fq_2_unpaired = genomeFileName + '_R2.s.fastq.gz'

    """
    trimmomatic \
    PE \
    -threads ${task.cpus} \
    -phred33 \
    ${genomeReads[0]} \
    ${genomeReads[1]} \
    $fq_1_paired \
    $fq_1_unpaired \
    $fq_2_paired \
    $fq_2_unpaired \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """

    stub:
    genomeFileName = genomeName.toString().split("\\_")[0]

    fq_1_paired = genomeFileName + '_R1.p.fastq.gz'
    fq_1_unpaired = genomeFileName + '_R1.s.fastq.gz'
    fq_2_paired = genomeFileName + '_R2.p.fastq.gz'
    fq_2_unpaired = genomeFileName + '_R2.s.fastq.gz'


    """
    echo "trimmomatic \
    PE \
    -threads ${task.cpus} \
    -phred33 \
    ${genomeReads[0]} \
    ${genomeReads[1]} \
    $fq_1_paired \
    $fq_1_unpaired \
    $fq_2_paired \
    $fq_2_unpaired \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36"
    
    touch ${genomeName}_R1.p.fastq.gz
    touch ${genomeName}_R2.p.fastq.gzcontainer 'quay.io/biocontainers/trimmomatic:0.35--6'
    """

}

workflow test {
    input_reads_ch = Channel.fromFilePairs("$launchDir/data/mock_data/*_{R1,R2}*fastq.gz")
    TRIMMOMATIC(input_reads_ch)
}
