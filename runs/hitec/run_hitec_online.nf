
params.data = '/data'
params.genomedir = '/data/genomes/'
//params.cfg = '-s 100m -k 23'
params.cfg = '-p 12'
// 'SRR065390'
orgTable = [
    'SaureusMW2'  : 'Staphylococcus aureus',
    'Hacinonychis':'Helicobacter acinonychis str. Sheeba'
]

genomeTable = [ 
    'SaureusMW2' :
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/265/GCF_000011265.1_ASM1126v1/GCF_000011265.1_ASM1126v1_genomic.fna.gz', 
    'Hacinonychis': //https://www.ncbi.nlm.nih.gov/genome/1444?genome_assembly_id=300646
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/305/GCF_000009305.1_ASM930v1/GCF_000009305.1_ASM930v1_genomic.fna.gz'
]

exptTable = [
    'SaureusMW2' : 'http://www.genomic.ch/edena/mw2Reads.seq.gz',
    'Hacinonychis':'http://sharcgs.molgen.mpg.de/data/reads_seq.gz'
]


paramsTable = [
    'SaureusMW2' : '  2820460  1.00 ',
    'Hacinonychis' : '  1557590   1.6 '
]

String[] parseExptID(String tx, String vx){
    tx.split(vx)
} 

orgIds = Channel.from(exptTable.keySet()).map{
    org -> [org, orgTable[org], genomeTable[org]]
}

process genomeDownload {
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

process bwamemIndex {
    tag{ orgDesc }

    storeDir "${params.genomedir}/bwa"

    input:
    set orgId, orgDesc, gnmFile from orgChan

    output:
    set orgId, orgDesc, gnmFile, file("${orgId}.fa.*") into idxChan

    """
    bwa index ${params.genomedir}/${orgId}.fa -p /tmp/${orgId}.fa
    cp /tmp/${orgId}.fa.bwt ${orgId}.fa.bwt
    cp /tmp/${orgId}.fa.pac ${orgId}.fa.pac
    cp /tmp/${orgId}.fa.ann ${orgId}.fa.ann
    cp /tmp/${orgId}.fa.bwt ${orgId}.fa.bwt
    cp /tmp/${orgId}.fa.amb ${orgId}.fa.amb
    cp /tmp/${orgId}.fa.sa ${orgId}.fa.sa
    rm /tmp/${orgId}.fa.*
    """

}

process SeqDownload {
    tag{ orgDesc }

    input:
    set orgId, orgDesc, gnmFile, idxFiles from idxChan

    output:
    set orgId, orgDesc, gnmFile, idxFiles, file("${exptId}.fasta") into pseqChan

    """
    wget ${exptTable[exptId]} -O ${exptId}.fasta.gz
    gunzip ${exptId}.fasta.gz
    """
}

(beforeChan1, beforeChan2) = pseqChan.into(2)

process runHiTEC{
    tag{ orgDesc }

    input:
    set orgId, orgDesc, gnmFile, idxFiles, file(beforeEC) from beforeChan1

    output:
    set orgId, orgDesc, gnmFile, idxFiles, file(beforeEC), file("afterEC.fasta") into ecChan

    //check params and root directory and fa/fq
    """
    hitec ${beforeEC} afterEC.fasta ${paramsTable[orgId]}
    """ 
}

process runBWABefore{
    tag{ orgDesc }

    input:
    set orgId, orgDesc, gnmFile, idxFiles, file(beforeEC) from beforeChan2

    output:
    set orgId, orgDesc, gnmFile, idxFiles, file(beforeEC), file("beforeEC.sam") into beforeSAMChan

    """
    bwa mem ${params.genomedir}/bwa/${orgId}.fa ${beforeEC} > beforeEC.sam
    """ 
}

process runBWAAfter{
    tag{ orgDesc }

    input:
    set orgId, orgDesc, gnmFile, idxFiles, file(beforeEC), file(afterEC) from ecChan

    output:
    set orgId, orgDesc, gnmFile, idxFiles, file(beforeEC), file("afterEC.sam") into afterSAMChan

    """
    bwa mem ${params.genomedir}/bwa/${orgId}.fa ${afterEC} > afterEC.sam
    """ 
}

mergedSAMChan = beforeSAMChan
    .join(afterSAMChan)
    .map {
        orgId1,
           orgDesc1, gnmFile1, idxFiles1, beforeEC, beforeSAM,
           orgDesc2, gnmFile2, idxFiles2, afterEC, afterSAM ->
        [ orgId1, orgDesc1, gnmFile1, idxFiles1,
          beforeEC, afterEC, beforeSAM, afterSAM  ]
    }

/*
processEvalEC{
    tag{ orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file(afterEC), file(beforeSAM), file(afterSAM) from mergedSAMChan

    output:
    set orgExptId, orgId, orgDesc, exptId, file("eval.json")

}
*/