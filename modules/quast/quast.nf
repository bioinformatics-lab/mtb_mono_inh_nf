
ch_refFILE = Channel.value("$baseDir/refFILE")

inputFilePattern = "./*_{R1,R2}.fastq.gz"
Channel.fromFilePairs(inputFilePattern)
        .into {  ch_in_PROCESS }



process process {
#    publishDir 'results/PROCESS'
#    container 'PROCESS_CONTAINER'


    input:
    set genomeFileName, file(genomeReads) from ch_in_PROCESS

    output:
    path("""${PROCESS_OUTPUT}""") into ch_out_PROCESS


    script:
    #FIXME
    genomeName= genomeFileName.toString().split("\\_")[0]
    
    """
    CLI PROCESS
    """
}



/*
#==============================================
# TODO quast
# quay.io/biocontainers/quast:5.0.2--1
#==============================================
*/

/*
   cp ../spades_results/work/ae/4d2677d866d9030d613d882de8b639/210_spades/scaffolds.fasta 210_scaffolds.fasta
   cp ../spades_results/work/b0/a4240ada1affe5b4217da78f21bc8f/23_spades/scaffolds.fasta 23_scaffolds.fasta
quast.py 23_scaffolds.fasta 210_scaffolds.fasta
*/
