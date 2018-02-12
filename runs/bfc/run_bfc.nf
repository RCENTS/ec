
params.data = '/data'
params.genomedir = '/data/genomes/'
//params.cfg = '-s 100m -k 23'
params.cfg = '-k 23'
// 'SRR065390'
org_table = [
    'E. coli K-12 MG1655': 'EcoliK12',
    'S. aureus MW2': 'SauresMW2'
]

genome_table = [ 
    'E. coli K-12 MG1655' :
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz', 
    'S. aureus MW2' :
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/265/GCF_000011265.1_ASM1126v1/GCF_000011265.1_ASM1126v1_genomic.fna.gz', 
]

// expt_table = [
//     'E. coli K-12 MG1655': ['SRR001665', 'SRR022918'],
//     'S. aureus MW2' : ['SRR022866']
// ]

expt_table = [
     'E. coli K-12 MG1655': ['SRR001665'],
    // 'S. aureus MW2' : ['SRR022918']
]

String[] parseExptID(String tx, String vx){
    tx.split(vx)
} 

org_ids = Channel.from(expt_table.keySet()).map{
    org -> [org, org_table[org], genome_table[org]]
}

process genomeDownload {
    tag{ org }

    storeDir params.genomedir

    input:
    set org, gnm, gurl from org_ids

    output:
    set org, gnm, file("${gnm}.fa") into org_chan

    """
    wget $gurl -O ${gnm}.fa.gz
    gunzip ${gnm}.fa.gz
    """
}

process bwamemIndex {
    tag{ org }

    storeDir "${params.genomedir}/bwa"

    input:
    set org, gnm, gfile from org_chan

    output:
    set org, gnm, gfile, file("${gnm}.fa.*") into idx_chan

    """
    bwa index ${params.genomedir}/${gnm}.fa -p /tmp/${gnm}.fa
    cp /tmp/${gnm}.fa.bwt ${gnm}.fa.bwt
    cp /tmp/${gnm}.fa.pac ${gnm}.fa.pac
    cp /tmp/${gnm}.fa.ann ${gnm}.fa.ann
    cp /tmp/${gnm}.fa.bwt ${gnm}.fa.bwt
    cp /tmp/${gnm}.fa.sa ${gnm}.fa.sa
    rm /tmp/${gnm}.fa.*
    """

}

// bwa_chan.subscribe{
//     println it
// }

expt_chan = idx_chan.flatMap {
    org, gnm, gnmfs, idxfs -> expt_table[org].collect {
        [org, gnm, gnmfs, idxfs, it]
    }
}.flatMap{
    org, gnm, gnmfs, idxfs, expt  ->
        sids = parseExptID(expt, ',')
        sids.collect{
            [org, gnm, gnmfs, idxfs, expt, it]
        }
}

process sraFetch {
    tag { org.toString() + " > " + expt.toString() + " > " + sid.toString() }

    input:
    set org, gnm, gnmfs, idxfs, expt, sid from expt_chan

    output:
    set org, gnm, gnmfs, idxfs, expt, sid, file("${sid}_*.fastq") into seq_chan

    """
    prefetch '$sid'
    vdb-validate '$sid'
    fastq-dump -I --split-files '$sid'
    """
}

process combinePairedEndFiles{
    tag { org.toString() + " > " + expt.toString() + " > " + sid.toString() }

    input:
    set org, gnm, gnmfs, idxfs, expt, sid, file(pseq) from seq_chan

    output:
    set org, gnm, gnmfs, idxfs, expt, sid, file("${sid}.fq") into fseq_chan
 
    """
    cat $pseq > ${sid}.fq
    """
}

oexp_chan = fseq_chan.map{ 
    org, gnm, gnmfs, idxfs, expt, sid, sfile -> 
        [org.toString() + "-" + expt.toString(),
         org, gnm, gnmfs, idxfs, expt, sid, sfile] 
}
.groupTuple()
.map{
    oeid, org, gnm, gnmfs, idxfs,  expt, sid, sfile -> 
        [oeid, org[0], gnm[0], gnmfs[0],
         idxfs[0], expt[0], sid, sfile]
}

// oexp_chan.subscribe{
//     println it
// }

process combineSRAs{
    tag { oeid.replace('-SRR', ' > SRR') }
    
    input:
    set oeid, org, gnm, gnmfs, idxfs, expt, sid, file(pseq) from oexp_chan

    output:
    set oeid, org, gnm, gnmfs, idxfs, expt, sid, file("before.fq") into cseq_chan

    script:
    if(sid.size() == 1)
        """
        mv $pseq before.fq
        """
    else
        """
        cat $pseq > before.fq
        """
}

process runBFC{
    tag { oeid.replace('-SRR', ' > SRR') }

    input:
    set oeid, org, gnm, gnmfs, idxfs, expt, sid, file(beforefs) from cseq_chan

    output:
    set oeid, org, gnm, gnmfs, idxfs, expt, sid, beforefs, file("after.fq") into ec_chan

    """
    bfc ${params.cfg} $beforefs > after.fq
    """ 
}


//full_seq.subscribe{ println "File: ${it.name}" }
