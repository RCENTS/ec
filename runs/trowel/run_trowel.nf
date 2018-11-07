

orgTable = [
    'EcoliK12MG1655'  : 'E. coli K-12 MG1655',
    'SauresMW2'       : 'S. aureus MW2',
    'ScerevisiaeS288C'   : 'S. cerevisiae S288C', 
    'Dmelanogaster' : 'D. melanogaster R6.18'
]

genomeTable = [ 
    'EcoliK12MG1655' :
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz', 
    'SauresMW2'      :
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/265/GCF_000011265.1_ASM1126v1/GCF_000011265.1_ASM1126v1_genomic.fna.gz', 
    'ScerevisiaeS288C'   :
      'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz',
    'Dmelanogaster' :
      'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/fasta/dmel-all-chromosome-r6.24.fasta.gz'
]

exptTable = [
    'EcoliK12MG1655' : [
                           'SRR001665', //D1
                           'ERR022075', //D2
                           'SRR022918'  //D3
                        ],
    'ScerevisiaeS288C' : ['SRR352384'], //D4
    'SauresMW2'      : ['SRR022866'], // D5
    'Dmelanogaster' : [
                           'SRR060098', //D6_1
                           'SRR018294', //D6_2
                           'SRR018292', //D6_3
                           'SRR018293'  //D6_4
                        ]
]

paramsTrowel = [ // parameters based on supplementary material
    'SRR001665' : '-k 19',
    'ERR022075' : '-k 19',
    'SRR022918' : '-k 19',
    'SRR352384' : '-k 19',
    'SRR022866' : '-k 17',
    'SRR060098' : '-k 19',
    'SRR018294' : '-k 19',
    'SRR018292' : '-k 19',
    'SRR018293' : '-k 19'
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

process GenomeIndexBWA {
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
    fastq-dump --split-files '${sraId}'
    """
}
//    fastq-dump -I --split-files '${sraId}'

(beforeChan1, beforeChan2) = pseqChan.into(2)

process Trowel{
    tag { orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC) from beforeChan1

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file("*.fastq.out") into ecChan

    """
    echo ${beforeEC} > inlist.txt
    trowel -f inlist.txt ${paramsTrowel[exptId]} -t ${params.nthreads} 
    """ 
}

process BWABeforeEC{
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

process BWAAfterEC{
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
    cp ${workflow.projectDir}/trowel.py .
    python trowel.py $beforeSAM $afterSAM ${orgExptId.replace(' ', '-')} > result1
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
    cp ${workflow.projectDir}/trowelbr.py .
    python trowelbr.py $beforeSAM $afterSAM ${orgExptId.replace(' ', '-')} > result2
    """

}

result_channel1.map{
    it.text
}.collectFile(name: 'trowel_eval_reads.tsv',
              storeDir: "${workflow.projectDir}",
              newLine: false)

result_channel2.map{
    it.text
}.collectFile(name: 'trowel_eval_bases.tsv',
              storeDir: "${workflow.projectDir}",
              newLine: false)

