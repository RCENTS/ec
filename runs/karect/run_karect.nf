
params.data = '/data'
params.genomedir = '/data/genomes'

//params.cfg = '-s 100m -k 23'

orgTable = [
    'h.pylori.454'  : 'h. pylori 454',
    'z.mobilis.454'     : 'z. mobilis 454',
    'c.elegans.ill'     : 'c. elegans ill'
    
]
genomeTable = [ 
    'h.pylori.454'  : 
        '',
    'z.mobilis.454'     : 
        '',
    'c.elegans.ill'     : 
        ' ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/sequence/genomic/c_elegans.PRJNA13758.WS241.genomic.fa.gz.'
]
exptTable = [
    'h.pylori.454'  : ['SRR023794', 'SRR023796'],
    'z.mobilis.454'     : ['SRR017972', 'SRR029606'],
    'c.elegans.ill'     : ['SRR065390 ']
]


String[] parseExptID(String tx, String vx){
    tx.split(vx)
} 

orgIds = Channel.from(exptTable.keySet()).map{
    org -> [org, orgTable[org], genomeTable[org]]
}
/*
params = Channel.from(exptTable.keySet()).map{
    par -> [org, paramsTessel[org], paramsBlue[org]]
}
*/

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

process bowtie2Index {
    tag{ orgDesc }

    storeDir "${params.genomedir}/bowtie2"

    input:
    set orgId, orgDesc, file(gnmFile) from orgChan

    output:
    set orgId, orgDesc, gnmFile, file("${orgId}.fa.*") into idxChan

    """
    bowtie2-build ${params.genomedir}/${orgId}.fa /tmp/${orgId}.fa
    cp /tmp/${orgId}.fa.1.bt2 ${orgId}.fa.1.bt2
    cp /tmp/${orgId}.fa.2.bt2 ${orgId}.fa.2.bt2
    cp /tmp/${orgId}.fa.3.bt2 ${orgId}.fa.3.bt2
    cp /tmp/${orgId}.fa.4.bt2 ${orgId}.fa.4.bt2
    cp /tmp/${orgId}.fa.rev.1.bt2 ${orgId}.fa.rev.1.bt2
    cp /tmp/${orgId}.fa.rev.2.bt2 ${orgId}.fa.rev.2.bt2
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
    set orgId, orgDesc, gnmFile, idxFiles, exptId, sraId, file("${sraId}.fastq") into fseqChan
 
    """
    cat ${pairedFiles} > ${sraId}.fastq
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
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file("beforeEC.fastq") into cseqChan

    script:
    if(sraIds.size() == 1)
        """
        mv ${sraFiles} beforeEC.fastq
        """
    else
        """
        cat ${sraFiles} > beforeEC.fastq
        """
}
(beforeChan1, beforeChan2) = cseqChan.into(2)

process runKarect{
    tag { orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC) from beforeChan1

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file("beforeEC_blue.fastq") into ecChan

    """
    ./karect -correct -threads=12 -matchtype=hamming -celltype=haploid -inputfile=${beforeEC}.fastq 
    """ 
}


process runBowtieBefore{
    tag { orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC) from beforeChan2

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file("beforeEC.sam") into beforeSAMChan

    """
    bowtie2 -x  ${params.genomedir}/bowtie2/${orgId}.fa  -U ${beforeEC} -S beforeEC.sam

    """ 
}

process runBowtieAfter{
    tag { orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file(afterEC) from ecChan

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file("afterEC.sam") into afterSAMChan

    """
    bowtie2 -x  ${params.genomedir}/bowtie2/${orgId}.fa -U ${afterEC} -S afterEC.sam
    """ 
}

mergedSAMChan = beforeSAMChan
    .merge(afterSAMChan)
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
