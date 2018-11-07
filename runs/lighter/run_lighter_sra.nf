

orgTable = [
    'EcoliK12MG1655'  : 'E. coli K-12 MG1655',
    'CelegansWS222'   : 'C. elegans WS222'
]

genomeTable = [ 
    'EcoliK12MG1655' :
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz', 
    'CelegansWS222'  :
      'ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/sequence/genomic/c_elegans.WS222.genomic.fa.gz'
]

exptTable = [
    'EcoliK12MG1655' : [
                           'ERR022075'
                        ],
     'CelegansWS222' : [
                           'SRR065390'
                        ]
]

paramsLighter = [ // parameters based on the paper
    'ERR022075' : '-K 19 4641650',
    'SRR065390' : '-K 23 100286000'
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

process GenomeIndexBowtie {
    tag{ orgDesc }

    storeDir "${params.genomedir}/bowtie2"

    input:
    set orgId, orgDesc, gnmFile from orgChan

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
        orgExptId = orgId.toString() + "-" + exptId.toString()
        sraIds.collect{
            [orgExptId, orgId, orgDesc, gnmFile, idxFiles , exptId, it]
        }
}
// .subscribe{
//     println it
// }

process FetchSRA {
    tag { orgId.toString() + " > " + exptId.toString() + " > " + sraId.toString() }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraId from exptChan

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, 
        exptId, sraId, file("${sraId}*.fastq") into pseqChan

    """
    prefetch ${sraId}
    vdb-validate '${sraId}'
    fastq-dump -I --split-files '${sraId}'
    """
}
//    fastq-dump -I --split-files '${sraId}'

(beforeChan1, beforeChan2) = pseqChan.into(2)

process Lighter{
    tag { orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC) from beforeChan1

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file("*.cor.*") into ecChan

    """
    lighter -r ${beforeEC} ${paramsLighter[exptId]} -t ${params.nthreads} 
    """ 
}

process BowtieBeforeEC{
    tag { orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC) from beforeChan2

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file("beforeEC.bam") into beforeSAMChan

    """
    bowtie2 -x  ${params.genomedir}/bowtie2/${orgId}.fa  -U ${beforeEC} | samtools view -bSh -F 0x900 - > bx.bam
    samtools sort -T bx.sorted -n -o beforeEC.bam bx.bam
    rm -rf bx.bam bx.sorted*
    """ 
}

process BowtieAfterEC{
    tag { orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file(afterEC) from ecChan

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(afterEC), file("afterEC.bam") into afterSAMChan

    """
    bowtie2 -x  ${params.genomedir}/bowtie2/${orgId}.fa  -U ${beforeEC} |  samtools view -bSh -F 0x900 - > ax.bam
    samtools sort -T ax.sorted -n -o afterEC.bam ax.bam
    rm -rf ax.bam ax.sorted*
    """ 
}

mergedSAMChan = beforeSAMChan
    .join(afterSAMChan)
    .map {
        orgExptId,
          orgId1, orgDesc1, gnmFile1, idxFiles1, exptId1, sraIds1, beforeEC, beforeSAM,
          orgId2, orgDesc2, gnmFile2, idxFiles2, exptId2, sraIds2, afterEC, afterSAM ->
        [ orgExptId, orgId1, orgDesc1, gnmFile1, idxFiles1, exptId1, sraIds1, 
          beforeEC, afterEC, beforeSAM, afterSAM  ]
    }

(mergedChan1, mergedChan2) = mergedSAMChan.into(2)

process EvalECReads{
    tag{ orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds,
        file(beforeEC), file(afterEC), file(beforeSAM), file(afterSAM) from mergedChan1

    output:
    file("result1") into result_channel1

    """
    cp ${workflow.projectDir}/lighterreads.py .
    python lighterreads.py $beforeSAM $afterSAM ${orgExptId.replace(' ', '-')} > result1
    """

}

process EvalECBases{
    tag{ orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds,
        file(beforeEC), file(afterEC), file(beforeSAM), file(afterSAM) from mergedChan2

    output:
    file("result2") into result_channel2

    """
    cp ${workflow.projectDir}/lighterbases.py .
    python lighterbases.py $beforeSAM $afterSAM ${orgExptId.replace(' ', '-')} > result2
    """

}

result_channel1.map{
    it.text
}.collectFile(name: 'lighter_eval_reads_sra.tsv',
              storeDir: "${workflow.projectDir}",
              newLine: false)

result_channel2.map{
    it.text
}.collectFile(name: 'lighter_eval_bases_sra.tsv',
              storeDir: "${workflow.projectDir}",
              newLine: false)

