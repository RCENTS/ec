orgTable = [
    'Llactis' : 'L. lactis', 
    'Tpallidum' : 'T. pallidum',
    'EColi' : 'E. coli K-12 MG1655',
    'Bsubtilis' : 'B. subtilis subtilis 168',
    'Saureus' : 'S. aureus MW2',
    'Paeruginosa' : 'P. aeruginosa PAO1',
    'Linterrogans4342' : 'L. interrogans str. 56601',
    'Linterrogans5823' : 'L. interrogans str. Fiocruz L1-130',
    'Hinfluenzae' : 'H. influenzae Rd KW20',
    'Scerevisiae' : 'S. cerevisiae S288C',
    'Celegans' : 'C. elegans WS222',
    'Dmelanogaster' : 'D. melanogaster R6.18'
]

genomeTable = [
    'Llactis' :
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/025/045/GCF_000025045.1_ASM2504v1/GCF_000025045.1_ASM2504v1_genomic.fna.gz',
    'Tpallidum' :
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/604/125/GCF_000604125.1_ASM60412v1/GCF_000604125.1_ASM60412v1_genomic.fna.gz',
    'EColi' : 
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz',
    'Bsubtilis' : 
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz',
    'Saureus' :
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/265/GCF_000011265.1_ASM1126v1/GCF_000011265.1_ASM1126v1_genomic.fna.gz',
    'Paeruginosa' :
      'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz',
    'Linterrogans4342' :
      'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/092/565/GCF_000092565.1_ASM9256v1/GCF_000092565.1_ASM9256v1_genomic.fna.gz',
    'Linterrogans5823' :
      'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/685/GCF_000007685.1_ASM768v1/GCF_000007685.1_ASM768v1_genomic.fna.gz',
    'Hinfluenzae' :
      'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/305/GCF_000027305.1_ASM2730v1/GCF_000027305.1_ASM2730v1_genomic.fna.gz',
    'Scerevisiae' : 
      'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz',
    'Celegans' :
      'ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/sequence/genomic/c_elegans.WS222.genomic.fa.gz',
    'Dmelanogaster' :
      'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/fasta/dmel-all-chromosome-r6.24.fasta.gz'
]

genomeNumTable = [
    'Llactis' : '2514220',
    'Tpallidum' : '1139180',
    'EColi' : '4641650',
    'Bsubtilis' : '4215610',
    'Saureus' : '2820460',
    'Paeruginosa' : '6264400',
    'Linterrogans4342' : '4698130',
    'Linterrogans5823' : '4627370',
    'Hinfluenzae' : '1830140',
    'Scerevisiae' : '12157100',
    'Celegans' : '102291899',
    'Dmelanogaster' : '120220296'
]

exptTable = [
    'Llactis' : ['SRR088759'],
    'Tpallidum' : ['SRR361468'],
    'EColi' : ['SRR001665', // SRR001665 === SRX000429
               'SRR396536', 
               'SRR396532',
               'SRR022918'],
    'Bsubtilis' : ['DRR000852'],
    'Saureus' : ['SRR022866'],
    'Paeruginosa' : ['SRR396641'],
    'Linterrogans4342' : ['SRR353563'],
    'Linterrogans5823' : ['SRR397962'],
    'Hinfluenzae' : ['SRR065202'],
    'Scerevisiae' : ['SRR352384'], // SRR352384 === SRX100885
    'Celegans' : ['SRR065390'],
    'Dmelanogaster' : [
        'SRR018292,SRR018293,SRR018294,SRR060098' 
            // SRR018292, SRR018293 === SRX006151 &&
            // SRR018294 === SRX006152 & 
            // SRR060098 === SRX023452
    ] 
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
//.subscribe{
//     println it
//}

process sraFetch {
    tag { orgDesc.toString() + " > " + exptId.toString() + " > " + sraId.toString() }

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
        [orgDesc.toString() + "-" + exptId.toString(),
         orgId, orgDesc, gnmFile, idxFiles, exptId, sraId, sraFile] 
}
.groupTuple()
.map{
    orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, sraFiles -> 
        [orgExptId, orgId[0], orgDesc[0], gnmFile[0],
         idxFiles[0], exptId[0], sraIds, sraFiles]
}

//oexptChan.subscribe{
//     println it
//}

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
    bwa mem ${params.genomedir}/bwa/${orgId}.fa ${beforeEC} | samtools view -bSh -F 0x900 - > bx.bam
    samtools sort -T bx.sorted -n -o beforeEC.bam bx.bam
    rm -rf bx.bam bx.sorted*
    """ 
}

process runBWAAfter{
    tag { orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file(afterEC) from ecChan

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(afterEC), file("afterEC.bam") into afterSAMChan

    """
    bwa mem ${params.genomedir}/bwa/${orgId}.fa ${afterEC} | samtools view -bSh -F 0x900 - > ax.bam
    samtools sort -T ax.sorted -n -o afterEC.bam ax.bam
    rm -rf ax.bam ax.sorted*
    """ 
}

mergedSAMChan = beforeSAMChan
    .join(afterSAMChan)
    .map {
        orgExptId,  // key
         orgId1, orgDesc1, gnmFile1, idxFiles1, exptId1, sraIds1, beforeEC, beforeSAM,
         orgId2, orgDesc2, gnmFile2, idxFiles2, exptId2, sraIds2, afterEC, afterSAM ->
        [ orgExptId, orgId1, orgDesc1, gnmFile1, idxFiles1,
          exptId1, sraIds1, beforeEC, afterEC, beforeSAM, afterSAM  ]
    }

// mergedSAMChan.subscribe{
//     println it
// }

(mergedChan1, mergedChan2) = mergedSAMChan.into(2)

process EvalEC{
    tag{ orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file(afterEC), file(beforeSAM), file(afterSAM) from mergedChan1

    output:
    file(result1) into result_channel1

    """
    cp ${workflow.projectDir}/racer.py .
    cp ${workflow.projectDir}/racerv2.py .
    python racerv2.py $beforeSAM $afterSAM ${orgExptId.replace(' ', '-')} > result1
    """ 

}

process readSearch{
    tag{ orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file(afterEC), file(beforeSAM), file(afterSAM) from mergedChan2

    output:
    file(result2) into result_channel2

    """
    printf ${orgExptId.replace(' ', '-')} > result2
    printf " " >> result2
    readSearch $gnmFile $beforeEC $afterEC | grep Gain | cut -d'=' -f2 | xargs >> result2
    """ 

}


result_channel1.map{
    it.text
}.collectFile(name: 'racer_eval.csv', 
              storeDir: "${workflow.projectDir}",
              newLine: false)


result_channel2.map{
    it.text
}.collectFile(name: 'racer_gain.txt', 
              storeDir: "${workflow.projectDir}",
              newLine: false)
