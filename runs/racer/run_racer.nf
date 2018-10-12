
params.data = '/project/schockalingam6/rcents/data'
params.genomedir = '/project/schockalingam6/rcents/data/genomes'

// 'SRR065390'
orgTable = [
    'Lactococcuslactis'  : 'Lactococcus lactis', 
    'Treponemapallidum'  : 'Treponema pallidum'
]

genomeTable = [ 
    'Lactococcuslactis' :
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/025/045/GCF_000025045.1_ASM2504v1/GCF_000025045.1_ASM2504v1_genomic.fna.gz',
    'Treponemapallidum' :
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/604/125/GCF_000604125.1_ASM60412v1/GCF_000604125.1_ASM60412v1_genomic.fna.gz'
]

exptTable = [
    'Lactococcuslactis' : ['SRR088759'],
    'Treponemapallidum' : ['SRR361468']

]
genomeNumTable = [
    'Lactococcuslactis'  : '2514220',
    'Treponemapallidum'  : '1139180'
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

process runRacer{
    tag { orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC) from beforeChan1

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file("afterEC.fastq") into ecChan

    """
    RACER ${beforeEC} afterEC.fastq ${genomeNumTable[orgId]}
    """ 
}

process runBWABefore{
    tag { orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC) from beforeChan2

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file("beforeEC.bam") into beforeSAMChan

    """
    bwa mem ${params.genomedir}/bwa/${orgId}.fa ${beforeEC} | samtools view -bS - > beforeEC.bam
    """ 
}

process runBWAAfter{
    tag { orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file(afterEC) from ecChan

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file("afterEC.bam") into afterSAMChan

    """
    bwa mem ${params.genomedir}/bwa/${orgId}.fa ${afterEC} | samtools view -bS - > afterEC.bam
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


process EvalEC{
    tag{ orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file(beforeSAM), file(afterSAM) from mergedSAMChan

    output:
    set file(result) into result_channel

    """
    wget https://raw.githubusercontent.com/RCENTS/ec/master/runs/racer/racer.py
    python racer.py $beforeSAM $afterSAM ${orgTable[orgId]} > result
    """ 

}

result_channel.map{
    it.text
}.collectFile(name: 'racer_results.txt', newLine: false)


