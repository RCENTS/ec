
// 'SRR065390'
orgTable = [
    'CelegansWS222'   : 'C. elegans WS222'
]

genomeTable = [ 
    'CelegansWS222'   :
      'ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/sequence/genomic/c_elegans.WS222.genomic.fa.gz'
]

exptTable = [
     'CelegansWS222' : ['SRR065390']
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

exptChan = idxChan.flatMap {
    orgId, orgDesc, gnmFile, idxFiles -> exptTable[orgId].collect {
        [orgId, orgDesc, gnmFile, idxFiles, it]
    }
}.flatMap{
    orgId, orgDesc, gnmFile, idxFiles , exptId  ->
        sraIds = parseExptID(exptId, ',')
        sraIds.collect{
            [orgId, orgDesc, gnmFile, idxFiles , exptId, it]
        }
}
// .subscribe{
//     println it
// }

process sraFetch {
    tag { orgId.toString() + " > " + exptId.toString() + " > " + sraId.toString() }

    input:
    set orgId, orgDesc, gnmFile, idxFiles, exptId, sraId from exptChan

    output:
    set orgId, orgDesc, gnmFile, idxFiles, exptId, sraId, file("${sraId}*.fastq") into pseqChan

    """
    prefetch ${sraId}
    vdb-validate '${sraId}'
    fastq-dump -I --split-files '${sraId}'
    """
}

process catPariedEndFiles{
    tag { orgId.toString() + " > " + exptId.toString() + " > " + sraId.toString() }

    input:
    set orgId, orgDesc, gnmFile, idxFiles, exptId, sraId, file(pairedFiles) from pseqChan

    output:
    set orgId, orgDesc, gnmFile, idxFiles, exptId, sraId, file("${sraId}.fq") into fseqChan
 
    """
    cat ${pairedFiles} > ${sraId}.fq
    """
}

oexptChan = fseqChan.map{ 
    orgId, orgDesc, gnmFile, idxFiles, exptId, sraId, sraFile -> 
        [orgId.toString() + "-" + exptId.toString(),
         orgId, orgDesc, gnmFile, idxFiles, exptId, sraId, sraFile] 
}
.groupTuple()
.map{
    orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, sraFiles -> 
        [orgExptId, orgId[0], orgDesc[0], gnmFile[0],
         idxFiles[0], exptId[0], sraIds, sraFiles]
}

// oexptChan.subscribe{
//     println it
// }

process catSRAFiles{
    tag { orgExptId.replace('-SRR', ' > SRR') }
    
    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(sraFiles) from oexptChan

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file("beforeEC.fq") into cseqChan

    script:
    if(sraIds.size() == 1)
        """
        mv ${sraFiles} beforeEC.fq
        """
    else
        """
        cat ${sraFiles} > beforeEC.fq
        """
}

(beforeChan1, beforeChan2) = cseqChan.into(2)

process runBFC{
    tag { orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC) from beforeChan1

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file("afterEC.fq") into ecChan

    """
    bfc ${params.cfg} ${beforeEC} > afterEC.fq
    """ 
}

process runBWABefore{
    tag { orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC) from beforeChan2

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file("beforeEC.sam") into beforeSAMChan

    """
    bwa mem ${params.genomedir}/bwa/${orgId}.fa ${beforeEC} > beforeEC.sam
    """ 
}

process runBWAAfter{
    tag { orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file(afterEC) from ecChan

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file("afterEC.sam") into afterSAMChan

    """
    bwa mem ${params.genomedir}/bwa/${orgId}.fa ${afterEC} > afterEC.sam
    """ 
}

mergedSAMChan = beforeSAMChan
    .join(afterSAMChan)
    .map {
        orgExptId1, orgId1, orgDesc1, gnmFile1, idxFiles1, exptId1, sraIds1, beforeEC1, beforeSAM,
        orgExptId2, orgId2, orgDesc2, gnmFile2, idxFiles2, exptId2, sraIds2, beforeEC2, afterSAM ->
        [ orgExptId1, orgId1, orgDesc1, gnmFile1, idxFiles1,
          exptId1, sraIds1, beforeEC1, beforeSAM, afterSAM  ]
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