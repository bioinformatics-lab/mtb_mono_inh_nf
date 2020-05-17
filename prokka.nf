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
