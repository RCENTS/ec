
params.data = '/data'
params.genomedir = '/data/genomes/'
//params.cfg = '-s 100m -k 23'
params.cfg = '-p 12'
// 'SRR065390'
orgTable = [
    'EcoliK12MG1655'  : 'E. coli K-12 MG1655',
]

genomeTable = [ 
    'EcoliK12MG1655' :
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz' 
]

exptTable = [
    'EcoliK12MG1655' : [
        'ERR008613', // 'ERA000206', 
        //'SRX00429',
        'SRR001665'
    ]
]

paramsTable = [
    'ERR008613' : ' 4641650   0.50 ',
    'SRR001665' : ' 4641650   0.44 '
]

String[] parseExptID(String tx, String vx){
    tx.split(vx)
} 

orgIds = Channel.from(exptTable.keySet()).map{
    org -> [org, orgTable[org], genomeTable[org]]
}

process GenomeDownload {
    tag{ orgDesc }

    storeDir params.genomedir

    input:
    set orgId, orgDesc, gnmURL from orgIds

    output:
    set orgId, orgDesc, file("${orgId}.fa") into orgChan

    """
    wget ${gnmURL} -O ${orgId}.fa.gz
    gunzip ${orgId}.fa.gz
    """
}

exptChan = orgChan.flatMap {
    orgId, orgDesc, gnmFile -> exptTable[orgId].collect {
        [orgId, orgDesc, gnmFile, it]
    }
}
// .subscribe{
//     println it
// }

process SRAFetch {
    tag { orgDesc.toString() + " > " + sraId.toString() }

    input:
    set orgId, orgDesc, gnmFile, sraId from exptChan

    output:
    set orgId, orgDesc, gnmFile, sraId, file("${sraId}.fastq") into pseqChan

    """
    prefetch ${sraId}
    vdb-validate '${sraId}'
    fastq-dump -I --split-files --gzip '${sraId}'
    gunzip -c ${sraId}_*.fastq.gz > ${sraId}.fastq
    rm -rf ${sraId}_1.fastq.gz
    rm -rf ${sraId}_2.fastq.gz
    """
}


//(beforeChan1, beforeChan2) = pseqChan.into(2)

process HiTEC{
    tag { orgDesc.toString() + " > " + sraId.toString() }

    input:
    set orgId, orgDesc, gnmFile, sraId, file(beforeEC) from pseqChan

    output:
    set orgId, orgDesc, gnmFile, sraId, file(beforeEC), file("afterEC.fastq") into ecChan

    //check params and root directory and fa/fq
    """
    hitec ${beforeEC} afterEC.fastq ${paramsTable[sraId]}
    """ 
}


process EvalEC{
    tag { orgDesc.toString() + " > " + sraId.toString() }

    input:
    set orgId, orgDesc, gnmFile, sraId, file(beforeEC), file(afterEC) from ecChan

    output:
    file(result) into result_channel

    """
    echo  "DATASET : " $sraId " ORGANISM : " $orgDesc > result
    readSearch $gnmFile $beforeEC $afterEC  >> result
    """ 

}


result_channel.map{
    it.text
}.collectFile(name: 'hitec_sra_gain.txt',
              storeDir: "${workflow.projectDir}",
              newLine: false)
