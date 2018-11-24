
orgTable = ['HSapiensChr14'   : 'H. sapiens chr14',
            'SAureusGAGE' : 'S. aureus GAGE']

dataTable = [
    'HSapiensChr14' : params.gagedataHS,
    'SAureusGAGE' : params.gagedataSA
]

genomeTable = [ 
    'HSapiensChr14'   :
      'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr14.fa.gz',
    'SAureusGAGE'  : // https://www.ncbi.nlm.nih.gov/genome/154?genome_assembly_id=299284
      'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/017/085/GCF_000017085.1_ASM1708v1/GCF_000017085.1_ASM1708v1_genomic.fna.gz'
]

paramsKarect = [
    'HSapiensChr14' : ' -matchtype=hamming -celltype=diploid ',
    'SAureusGAGE' : ' -matchtype=hamming -celltype=haploid '
]

paramsKarectAlign  = [
    'HSapiensChr14' : ' -matchtype=hamming  ',
    'SAureusGAGE' : ' -matchtype=hamming  '
]

orgIds = Channel.from(orgTable.keySet()).map{
    org -> [org, orgTable[org], genomeTable[org] ]
}

process GenomeDownload {
    tag{ orgDesc }

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

(beforeChan1, beforeChan2) = orgChan.into(2)


process KarectAlignBefore{
    tag { orgDesc }

    input:
    set orgId, orgDesc, gnmFile from beforeChan1

    output:
    set orgId, orgDesc, gnmFile, file("beforeEC.txt") into beforeSAMChan

    """
    karect -align -threads=${params.nthreads}  ${paramsKarectAlign[orgId]}  -refgenomefile=$gnmFile  -inputfile=${params.gagedatadir}/${dataTable[orgId]} -alignfile=beforeEC.txt
    """ 

}

process Karect{
    tag { orgDesc }

    input:
    set orgId, orgDesc, gnmFile from beforeChan2

    output:
    set orgId, orgDesc, gnmFile, file("afterEC.fastq") into ecChan


    """
    karect -correct -threads=${params.nthreads} ${paramsKarect[orgId]} -inputfile=${beforeEC} 
    """ 
}


mergedSAMChan = beforeSAMChan
    .join(ecChan)
    .map {
        orgId1, 
            orgDesc1, gnmFile1, beforeSAM,
            orgDesc2, gnmFile2, afterEC ->
        [ orgId1, orgDesc1, gnmFile1, idxFiles1, beforeSAM, afterEC ]
    }



process EvalEC{
    tag{ orgDesc }

    input:
    set orgId, orgDesc, gnmFile, 
        file(beforeEC),  file(beforeSAM), file(afterEC) from mergedSAMChan

    output:
    file("result1") into result_channel1


    """
    karect -eval -threads=${params.nthreads}  ${paramsKarectAlign[orgId]}  -inputfile=${params.gagedatadir}/${dataTable[orgId]}  -resultfile=$afterEC -refgenomefile=$gnmFile -alignfile=$beforeSAM -evalfile=result.txt
    echo "DATASET : GAGE   ORGANISM : " $orgDesc  > result1
    tail -n 20 result.txt >> result1
    rm result.txt
    """

}

result_channel1.map{
    it.text
}.collectFile(name: 'karect_eval_gage.txt',
              storeDir: "${workflow.projectDir}",
              newLine: false)
