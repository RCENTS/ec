
params.data = '/data'
params.cfg = '-t8 -s 100m -k 23'
// 'SRR065390'
sra_ids = Channel.from('SRR001665', 'SRR022918')
// env {
//     NXF_WORK = '/data'
// }
process sraFetch {
    input:
    val sid from sra_ids

    output:
    set val("$sid"), file("${sid}_*.fastq") into sequences

    """
    prefetch '$sid'
    vdb-validate '$sid'
    fastq-dump -I --split-files '$sid'
    """

}

    // """
    // echo 1 > ${sid}_1.fastq
    // echo 2 > ${sid}_2.fastq
    // """

process combineFiles{
    input:
    set val(sid), file(seq) from sequences

    output:
    set val("$sid"), file("${sid}.fq") into full_seq

    """
    cat $seq > ${sid}.fq
    """

}
//sequences
   // .flatMap()
//    .subscribe { println "${it}" }

process runBFC{
    input:
    set val(sid), file(fseq) from full_seq

    output:
    set val("$sid"), file('bfc_out') into bfc_rout

    """
    bfc ${params.cfg} $fseq > bfc_out
    """ 
}
//full_seq.subscribe{ println "File: ${it.name}" }