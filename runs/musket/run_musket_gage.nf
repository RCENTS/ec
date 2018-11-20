
orgTable = [
    'HSapiensChr14'   : 'H. sapiens (chr14)'
]

genomeTable = [ 
    'HSapiensChr14'   :
      'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr14.fa.gz'
]


orgIds = Channel.from(orgTable.keySet()).map{
    org -> [org, orgTable[org], genomeTable[org] ]
}

process GenomeDownload {
    tag{ 'H. sapiens (chr14)' }

    storeDir params.genomedir

    input:
    set orgId, orgDesc, gnmURL  from orgIds

    output:
    set orgId, orgDesc, file("${orgId}.fa") into orgChan

    """
    wget ${gnmURL} -O ${orgId}.fa.gz
    gunzip ${orgId}.fa.gz
    """
}

process BWAIndex {
    tag{ 'H. sapiens (chr14)' }

    storeDir "${params.genomedir}/bwa"

    input:
    set orgId, orgDesc, gnmFile from orgChan

    output:
    set orgId, orgDesc,  gnmFile, file("${orgId}.fa.*") into idxChan

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

(beforeChan1, beforeChan2) = idxChan.into(2)


process BWABefore{
    tag { 'BWA Before Musket' }

    input:
    set orgId, orgDesc, gnmFile, idxFiles from beforeChan1

    output:
    set orgId, orgDesc, gnmFile, idxFiles, file("beforeEC.fastq"), file("beforeEC.sam") into beforeSAMChan

    """
    cat ${params.gagedatadir}/${params.gagedata} | paste - - - - | sort -k1,1 -t " " -S 64G -T ./ --parallel=${sortThreads} | tr "\t" "\n" > beforeEC.fastq
    bwa aln -t ${params.nthreads} -n ${params.hammingdist} -o 0 ${params.genomedir}/bwa/${orgId}.fa ${params.gagedatadir}/${params.gagedata}  > beforeEC.sai
    bwa samse ${params.genomedir}/bwa/${orgId}.fa  beforeEC.sai  ${beforeEC} | samtools view -bSh -F 0x900 - > beforeEC.bam
    samtools view -H beforeEC.bam > beforeEC.sam
    samtools view beforeEC.bam | sort -k1,1 -t ' ' -T ./ -S64G  >>  beforeEC.sam
    """ 

}

process Musket{
    tag { 'Musket Run' }

    input:
    set orgId, orgDesc, gnmFile, idxFiles from beforeChan2

    output:
    set orgId, orgDesc, gnmFile, idxFiles, file("afterEC.fastq") into ecChan


   script:
    if (sortThreads < 2)
        """
        musket -p ${params.nthreads} ${params.gagedatadir}/${params.gagedata} -o tmp.fastq
        cat tmp.fastq | paste - - - - | sort -k1,1 -t " " -S 64G -T ./ | tr "\t" "\n" > afterEC.fastq
        rm tmp.fastq
        """ 
    else
        """
        musket -p ${params.nthreads}  ${params.gagedatadir}/${params.gagedata} -o tmp.fastq
        cat tmp.fastq | paste - - - - | sort -k1,1 -t " " -S 64G -T ./ --parallel=${sortThreads} | tr "\t" "\n" > afterEC.fastq
        rm tmp.fastq
        """ 
}


mergedSAMChan = beforeSAMChan
    .join(ecChan)
    .map {
        orgId1, 
            orgDesc1, gnmFile1, idxFiles1, beforeEC, beforeSAM,
            orgDesc2, gnmFile2, idxFiles2, afterEC ->
        [ orgId1, orgDesc1, gnmFile1, idxFiles1, beforeEC, beforeSAM, afterEC ]
    }



process EvalEC{
    tag{ 'Musket Eval Read' }

    input:
    set orgId, orgDesc, gnmFile, idxFiles, 
        file(beforeEC),  file(beforeSAM), file(afterEC) from mergedChan1

    output:
    file("result1") into result_channel1

    """
    sam-analysis.py -a ambig.lst -t alntrim.lst $beforeSAM align.tef unmapped.lst
    quake-analy.py -f $beforeEC -c $afterEC -o correction.tef -t cortrim.lst
    cat ambig.lst unmapped.lst | sort > ua.lst
    comp2pcalign correction.tef align.tef ua.lst ${params.hammingdist} result.txt alntrim.lst 
    echo "DATASET : " $exptId  " ORGANISM : " $orgDesc  > result1
    cat result.txt >> result1
    """

}

result_channel1.map{
    it.text
}.collectFile(name: 'musket_eval_gage.txt',
              storeDir: "${workflow.projectDir}",
              newLine: false)
