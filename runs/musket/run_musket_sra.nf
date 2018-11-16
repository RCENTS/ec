
orgTable = [
    'CelegansWS222'   : 'C. elegans',
    'EcoliK12MG1655'  : 'E. coli K-12 MG1655'
] 

genomeTable = [ 
    'CelegansWS222' :
      'ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/sequence/genomic/c_elegans.WS222.genomic.fa.gz',
    'EcoliK12MG1655' :
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz' 
]

exptTable = [
//    'CelegansWS222'  : ['SRR065390'],
    'EcoliK12MG1655' : ['ERR022075']
]

String[] parseExptID(String tx, String vx){
    tx.split(vx)
} 

orgIds = Channel.from(exptTable.keySet()).map{
    org -> [org, orgTable[org], genomeTable[org]]
}

sortThreads = params.nthreads - 3 

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

process BWAIndex {
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

process SRAFetch {
    tag { orgId.toString() + " > " + exptId.toString() + " > " + sraId.toString() }

    input:
    set orgId, orgDesc, gnmFile, idxFiles, exptId, sraId from exptChan

    output:
    set orgId, orgDesc, gnmFile, idxFiles, exptId, sraId, file("${sraId}*.fastq.gz") into pseqChan

    """
    prefetch ${sraId}
    vdb-validate '${sraId}'
    fastq-dump -I --split-files --gzip '${sraId}'
    """
}

process ConcatPariedEndFiles{
    tag { orgId.toString() + " > " + exptId.toString() + " > " + sraId.toString() }

    input:
    set orgId, orgDesc, gnmFile, idxFiles, exptId, sraId, file(pairedFiles) from pseqChan

    output:
    set orgId, orgDesc, gnmFile, idxFiles, exptId, sraId, file("${sraId}.fastq") into fseqChan

    script: 
    if (sortThreads < 2)
        """
        gunzip -c ${pairedFiles} | paste - - - - | sort -k1,1 -t ' ' -S 64G | tr '\t' '\n' > ${sraId}.fastq
        """ 
    else
        """
        gunzip -c ${pairedFiles} | paste - - - - | sort -k1,1 -t ' ' -S 64G --parallel=${sortThreads} | tr '\t' '\n' > ${sraId}.fastq
        """ 
}

oexptChan = fseqChan.map{ 
    orgId, orgDesc, gnmFile, idxFiles, exptId, sraId, sraFile -> 
        [orgId.toString() + "-" + exptId.toString(),
         orgId, orgDesc, gnmFile, idxFiles, exptId, sraId, sraFile] 
}


(beforeChan1, beforeChan2) = oexptChan.into(2)

process Musket{
    tag { orgId.toString() + " > " + exptId.toString() }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC) from beforeChan1

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file("afterEC.fastq") into ecChan

    script:
    if (sortThreads < 2)
        """
        musket -p ${params.nthreads}  ${beforeEC} -o tmp.fastq
        cat tmp.fastq | paste - - - - | sort -k1,1 -t " " -S 64G | tr "\t" "\n" > afterEC.fastq
        rm tmp.fastq
        """ 
    else
        """
        musket -p ${params.nthreads}  ${beforeEC} -o tmp.fastq
        cat tmp.fastq | paste - - - - | sort -k1,1 -t " " -S 64G --parallel=${sortThreads} | tr "\t" "\n" > afterEC.fastq
        rm tmp.fastq
        """ 
}


process BWABefore{
    tag { orgId.toString() + " > " + exptId.toString() }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC) from beforeChan2

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file("beforeEC.sam") into beforeSAMChan

    """
    bwa aln -t ${params.nthreads} -n 4 -o 0 ${params.genomedir}/bwa/${orgId}.fa ${beforeEC} > beforeEC.sai
    bwa samse -t ${params.nthreads}  ${params.genomedir}/bwa/${orgId}.fa  beforeEC.sai | samtools view -bSh -F 0x900 - > bx.bam
    samtools sort -T bx.sorted -n -o beforeEC.bam bx.bam
    rm -rf bx.bam bx.sorted* beforeEC.sai 
    samtools view -h  -o beforeEC.sam
    """ 
}

mergedSAMChan = beforeSAMChan
    .join(ecChan)
    .map {
        orgExptId,
            orgId1, orgDesc1, gnmFile1, idxFiles1, exptId1, sraIds1, beforeEC, beforeSAM,
            orgId2, orgDesc2, gnmFile2, idxFiles2, exptId2, sraIds2, afterEC ->
        [ orgExptId,
            orgId1, orgDesc1, gnmFile1, idxFiles1, exptId1, sraIds1,
            beforeEC, beforeSAM, afterEC  ]
    }

process EvalEC{
    tag{ orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds,
        file(beforeEC), file(beforeSAM), file(afterEC) from mergedSAMChan

    output:
    file("result.txt") into result_channel1

    """
    cat x > result.txt
    """

}
/**/
