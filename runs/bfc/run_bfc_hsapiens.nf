
// 'SRR065390'
orgTable = [
    'HSapiens37d5'   : 'H. sapiens hs37d5'
]

genomeTable = [ 
    'HSapiens37d5'   :
      'ftp://ftp.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz'
]


exptTable = [
     'HSapiens37d5' : [params.hs37d5data]
]


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
    orgId, orgDesc, gnmFile, idxFiles -> 
        exptFile = exptTable[orgId]
        exptFilePaths = parseExptID(exptFile, ',')
        [orgId, orgDesc, gnmFile, idxFiles , exptFile, exptFilePaths]
    
}

process catExptFiles{
    tag { 'H Sapiens 37d5 cat' }

    storeDir params.hs37d5datadir

    input:
    set orgId, orgDesc, gnmFile, idxFiles, exptFile, file(exptFilePaths) from oexptChan

    output:
    set orgId, orgDesc, gnmFile, idxFiles, exptFile, file("HSapiens37d5.fastq.gz") into cseqChan

    script:
    if(exptFilePaths.size() == 1)
        """
        mv ${exptFilePaths} HSapiens37d5.fastq.gz
        """
    else
        """
        gunzip -c ${exptFilePaths} | gzip -1 > HSapiens37d5.fastq.gz
        """
}

(beforeChan1, beforeChan2) = cseqChan.into(2)


process runBFC{
    tag { 'BFC Run' }

    storeDir params.hs37d5datadir

    input:
    set orgId, orgDesc, gnmFile, idxFiles, file(beforeEC) from beforeChan1

    output:
    set orgId, orgDesc, gnmFile, idxFiles, file("HSapiens37d5.fastq.gz") into ecChan

    """
    bfc ${params.cfg} ${beforeEC} | gzip -1 > HSapiens37d5-BFC.fastq.gz
    """ 
}

 process randomSelect {
    tag { 'Random Select Before BFC' }

    input:
    set orgId, orgDesc, gnmFile, idxFiles, file(beforeEC) from beforeChan2

    output:
    set orgId, orgDesc, gnmFile, idxFiles, file(beforeEC), file("beforeRSelect.fastq.gz"), file("beforeRSIds.txt") into rselChan

    """
    seqtk sample -s100 ${beforeEC} ${params.nsample} | gzip -1 > beforeRSelect.fastq.gz
    gunzip -c beforeRSelect.fastq.gz | awk '/^@/ {print $1}' | cut -d'@' -f2 > beforeRSIds.txt
    """
}

process runBWABefore{
    tag { 'BWA Before BFC' }

    input:
    set orgId, orgDesc, gnmFile, idxFiles, file(beforeEC), file(beforeRSel) file(beforeRSelIds) from rselChan

    output:
    set orgId, orgDesc, gnmFile, idxFiles, file(beforeEC), file(beforeRSel) file(beforeRSelIds), file("beforeEC.bam") into beforeSAMChan

    """
    bwa mem ${params.genomedir}/bwa/${orgId}.fa ${beforeEC} | samtools view -bSh -F 0x900 - > bx.bam
    samtools sort -T bx.sorted -n -o beforeEC.bam bx.bam 
    rm -f bx.bam bx.sorted*
    """ 
}

mergedSAMChan = beforeSAMChan
    .join(ecChan)
    .map {
        orgId1, 
            orgDesc1, gnmFile1, idxFiles1, beforeEC, beforeRSel, beforeRSelIds, beforeSAM,
            orgDesc2, gnmFile2, idxFiles2, afterEC ->
        [ orgId1, orgDesc1, gnmFile1, idxFiles1,
          beforeEC, beforeRSel, beforeRSelIds, beforeSAM, afterEC  ]
    }


process runBWAAfter{
    tag { 'BWA After BFC' }

    input:
    set orgId, orgDesc, gnmFile, idxFiles, beforeEC, 
        file(beforeRSel), file(beforeRSelIds), file(beforeSAM), file(afterEC) from mergedSAMChan

    output:
    set orgId, orgDesc, gnmFile, idxFiles, beforeEC,
        beforeRSel, beforeRSelIds, beforeSAM, 
        file(afterEC), file("afterRSelect.fastq.gz"), file("afterEC.bam") into afterSAMChan

    """
    seqtk subseq ${afterEC} ${beforeRSelIds} | gunzip -1 > afterRSelect.fastq.gz
    bwa mem ${params.genomedir}/bwa/${orgId}.fa afterRSelect.fastq.gz| samtools view -bSh -F 0x900 - > ax.bam
    samtools sort -T ax.sorted -n -o afterEC.bam ax.bam
    rm -f ax.bam ax.sorted*
    """ 
}


process EvalEC{
    tag { 'H Sapies 37d5 eval' }

    input:
    set orgId, orgDesc, gnmFile, idxFiles, beforeEC,
        beforeRSel, beforeRSelIds, file(beforeSAM),
        file(afterEC), file(afterRSel), file(afterSAM) from mergedSAMChan

    output:
    file("result1") into result_channel1

    """
    cp ${workflow.projectDir}/bfc.py .
    python bfc.py $beforeSAM $afterSAM ${orgExptId.replace(' ', '-')} > result1
    """ 

}

result_channel1.map{
    it.text
}.collectFile(name: 'bfc_eval_hsapiens.tsv', 
              storeDir: "${workflow.projectDir}",
              newLine: false)

