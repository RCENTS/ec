
orgTable = [
    'HSapiensChr14'   : 'H. sapiens (chr14)'
]

genomeTable = [ 
    'HSapiensChr14'   :
      'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr14.fa.gz'
]


exptTable = [
     'HSapiensChr14' : params.gagedata
]


orgIds = Channel.from(exptTable.keySet()).map{
    org -> [org, orgTable[org], genomeTable[org], exptTable[org]]
}

process GenomeDownload {
    tag{ 'H. sapiens (chr14)' }

    storeDir params.genomedir

    input:
    set orgId, orgDesc, gnmURL, exptFile from orgIds

    output:
    set orgId, orgDesc, exptFile, file("${orgId}.fa") into orgChan

    """
    wget ${gnmURL} -O ${orgId}.fa.gz
    gunzip ${orgId}.fa.gz
    """
}

process GenomeIndexBowtie {
    tag{ 'H. sapiens (chr14)' }

    storeDir "${params.genomedir}/bowtie2"

    input:
    set orgId, orgDesc, exptFile, gnmFile from orgChan

    output:
    set orgId, orgDesc, exptFile, gnmFile, file("${orgId}.fa.*") into idxChan

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

(beforeChan1, beforeChan2) = idxChan.into(2)


process BowtieBefore{
    tag { 'BWA Before Lighter' }

    input:
    set orgId, orgDesc, file(beforeEC), gnmFile, idxFiles from beforeChan1

    output:
    set orgId, orgDesc, gnmFile, idxFiles, file(beforeEC), file("beforeEC.bam") into beforeSAMChan

    """
    bowtie2 -x  ${params.genomedir}/bowtie2/${orgId}.fa  -U ${beforeEC} | samtools view -bSh -F 0x900 - > bx.bam
    samtools sort -T bx.sorted -n -o beforeEC.bam bx.bam
    rm -rf bx.bam bx.sorted*
    """ 
}

process Lighter{
    tag { 'Lighter Run' }

    storeDir params.gagedatadir

    input:
    set orgId, orgDesc, file(beforeEC), gnmFile, idxFiles from beforeChan1

    output:
    set orgId, orgDesc, gnmFile, idxFiles, file("*.cor.*") into ecChan

    """
    lighter -r ${beforeEC} -K 19  107043718 -t ${params.nthreads}
    """ 
}

process BowtieAfter{
    tag { 'Bowtie After Lighter' }

    input:
    set orgId, orgDesc, gnmFile, idxFiles, file(afterEC) from ecChan

    output:
    set orgId, orgDesc, gnmFile, idxFiles, file(afterEC), file("afterEC.bam") into afterSAMChan

    """
    bowtie2 -x  ${params.genomedir}/bowtie2/${orgId}.fa  -U ${afterEC} |  samtools view -bSh -F 0x900 - > ax.bam
    samtools sort -T ax.sorted -n -o afterEC.bam ax.bam
    rm -rf ax.bam ax.sorted*
    """ 
}


mergedSAMChan = beforeSAMChan
    .join(afterSAMChan)
    .map {
        orgId1, 
            orgDesc1, gnmFile1, idxFiles1, beforeEC, beforeSAM,
            orgDesc2, gnmFile2, idxFiles2, afterEC, afterSAM ->
        [ orgId1, orgDesc1, gnmFile1, idxFiles1,
          beforeEC, beforeSAM, afterEC, afterSAM  ]
    }


(mergedChan1, mergedChan2) = mergedSAMChan.into(2)


process EvalECReads{
    tag{ 'Lighter Eval Read' }

    input:
    set orgId, orgDesc, gnmFile, idxFiles, 
        file(beforeEC), file(beforeSAM), file(afterEC), file(afterSAM) from mergedChan1

    output:
    file("result1") into result_channel1

    """
    cp ${workflow.projectDir}/lighterreads.py .
    python lighterreads.py $beforeSAM $afterSAM GAGE > result1
    """

}

process EvalECBases{
    tag{ 'Lighter Eval Bases' }

    input:
    set orgId, orgDesc, gnmFile, idxFiles, 
        file(beforeEC), file(beforeSAM), file(afterEC), file(afterSAM) from mergedChan2

    output:
    file("result2") into result_channel2

    """
    cp ${workflow.projectDir}/lighterbases.py .
    python lighterbases.py $beforeSAM $afterSAM GAGE > result2
    """

}

result_channel1.map{
    it.text
}.collectFile(name: 'lighter_eval_reads_gage.tsv',
              storeDir: "${workflow.projectDir}",
              newLine: false)

result_channel2.map{
    it.text
}.collectFile(name: 'lighter_eval_bases_gage.tsv',
              storeDir: "${workflow.projectDir}",
              newLine: false)
