
/*
unzip dados_emilyn-20200315T045112Z-001.zip  -d ./001/
unzip dados_emilyn-20200315T045112Z-002.zip  -d ./002/
unzip dados_emilyn-20200315T045112Z-003.zip  -d ./003/
unzip dados_emilyn-20200315T045112Z-004.zip  -d ./004/

rm *zip

mv 001/dados_emilyn/* .
mv 002/dados_emilyn/* .
mv 003/dados_emilyn/* .
mv 004/dados_emilyn/* .

rm -rf 001 002 003 004

mv 2903_S1_L001_R1_001.fastq.gz 2703_S1_L001_R1_001.fastq.gz
mv 2903_S1_L001_R2_001.fastq.gz 2703_S1_L001_R2_001.fastq.gz


mv 2703_S1_L001_R1_001.fastq.gz 2703_R1.fastq.gz
mv 2703_S1_L001_R2_001.fastq.gz 2703_R2.fastq.gz
mv 3384_S2_L001_R1_001.fastq.gz 3384_R1.fastq.gz
mv 3384_S2_L001_R2_001.fastq.gz 3384_R2.fastq.gz
mv 3652_S7_L001_R1_001.fastq.gz 3652_R1.fastq.gz
mv 3652_S7_L001_R2_001.fastq.gz 3652_R2.fastq.gz
mv 4025_S5_L001_R1_001.fastq.gz 4025_R1.fastq.gz
mv 4025_S5_L001_R2_001.fastq.gz 4025_R2.fastq.gz
mv 4084_S4_L001_R1_001.fastq.gz 4084_R1.fastq.gz
mv 4084_S4_L001_R2_001.fastq.gz 4084_R2.fastq.gz
mv 4106_S9_L001_R1_001.fastq.gz 4106_R1.fastq.gz
mv 4106_S9_L001_R2_001.fastq.gz 4106_R2.fastq.gz
mv 4736_S10_L001_R1_001.fastq.gz 4736_R1.fastq.gz
mv 4736_S10_L001_R2_001.fastq.gz 4736_R2.fastq.gz
mv 5765_S3_L001_R1_001.fastq.gz 5765_R1.fastq.gz
mv 5765_S3_L001_R2_001.fastq.gz 5765_R2.fastq.gz
mv 5800_S8_L001_R1_001.fastq.gz 5800_R1.fastq.gz
mv 5800_S8_L001_R2_001.fastq.gz 5800_R2.fastq.gz
mv 5923_S6_L001_R1_001.fastq.gz 5923_R1.fastq.gz
mv 5923_S6_L001_R2_001.fastq.gz 5923_R2.fastq.gz

*/

/*
################
NEXTFLOW Global Config
################
*/

params.outdir = "results"
ch_refGbk = Channel.value("$baseDir/NC000962_3.gbk")
ch_refFasta = Channel.value("$baseDir/NC000962_3.fasta")

Channel.fromFilePairs("./*_{R1,R2}.fastq.gz", flat: true)
        .into { ch_fastqGz; ch_snippy; ch_tbProfiler }


/*
################
gzip these files
################
*/


process gzip {
    container 'centos:8'

    input:
    tuple (genomeName, path(read_1_gz), path(read_2_gz)) from ch_fastqGz

    output:
    tuple genomeName, path(genome_1_fq), path(genome_2_fq) into ch_trimmomatic
    tuple genomeName, path(genome_1_fq), path(genome_2_fq) into ch_in_rdAnalyzer
    tuple genomeName, path(genome_1_fq), path(genome_2_fq) into ch_in_spades

    script:
    genome_1_fq = read_1_gz.name.split("\\.")[0] + '.fastq'
    genome_2_fq = read_2_gz.name.split("\\.")[0] + '.fastq'
    """
    gzip -dc -k ${read_1_gz} > ${genome_1_fq}
    gzip -dc -k ${read_2_gz} > ${genome_2_fq}
    """

}


/*
###############
trimmomatic
###############
*/


process trimmomatic {
    container 'quay.io/biocontainers/trimmomatic:0.35--6'


    input:
    tuple genomeName, path(fq_1), path(fq_2) from ch_trimmomatic

    output:
    tuple genomeName, path(fq_1_paired) into ch_in_spotyping

    script:
    fq_1_paired = genomeName + '_1_paired.fastq'
    fq_1_unpaired = genomeName + '_1_unpaired.fastq' //single
    fq_2_paired = genomeName + '_2_paired.fastq'
    fq_2_unpaired = genomeName + '_2_unpaired.fastq' //single


    """
    trimmomatic \
    PE -phred33 \
    $fq_1 \
    $fq_2 \
    $fq_1_paired \
    $fq_1_unpaired \
    $fq_2_paired \
    $fq_2_unpaired \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """
}

/*
#==============================================
# TB_Profiler
#==============================================
*/


process tbProfiler {

    container 'quay.io/biocontainers/tb-profiler:2.8.6--pypy_0'

    input:
    tuple (genomeName, path(read_1_gz), path(read_2_gz)) from ch_tbProfiler

    script:

    """
    tb-profiler profile -1 $read_1_gz -2 $read_2_gz  -t 4 -p $genomeName
    """
}


/*
###############
Spotyping
###############
*/



process spotyping {
    container 'abhi18av/spotyping'

    input:
    tuple genomeName, path(fq_1_paired) from ch_in_spotyping

    script:

    """
    python /SpoTyping-v2.0/SpoTyping-v2.0-commandLine/SpoTyping.py ./${fq_1_paired} -o ${genomeName}.txt
    """
}


///*
//###############
//Spades
//- Run on cloud
//- Extremely RAM heavy
//###############
//*/
//
//process spades {
//    container 'quay.io/biocontainers/spades:3.14.0--h2d02072_0'
//
//    input:
//    tuple genomeName, path(fq_1), path(fq_2) from ch_in_spades
//
//    script:
//
//
//    """
//    spades.py -k 21,33,55,77 --careful --only-assembler --pe1-1 ${fq_1} --pe1-2 ${fq_2} -o ${genomeName}_spades -t 2
//
//    """
//}


#!/usr/bin/env nextflow


/*
#==============================================
# config
#==============================================
*/

params.outdir = "results"
ch_refGbk = Channel.value("$baseDir/NC000962_3.gbk")
ch_refFasta = Channel.value("$baseDir/NC000962_3.fasta")


/*
#==============================================
# read genomes
#==============================================
*/

// TODO update this to receive the output of spades process
Channel.fromPath("./*_scaffolds.fasta")
        .into { ch_in_prokka }


/*
#==============================================
# prokka
#==============================================
*/


process prokka {
   container 'quay.io/biocontainers/prokka:1.14.6--pl526_0'
   publishDir 'results/prokka'

echo true

   input:
   path refFasta from ch_refFasta
   path bestContig from ch_in_prokka

   script:
   genomeName = bestContig.getName().split("\\_")[0]

   """
   prokka --outdir ./${genomeName}_prokka --prefix $genomeName ${bestContig}
   """

}

