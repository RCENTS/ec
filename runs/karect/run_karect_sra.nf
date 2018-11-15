
orgTable = [
    'HPylori'  : 'H. pylori v225d',
    'Zmobilis'     : 'Z. mobilis subsp. mobilis ATCC 31821',
    'CelegansWS241'     : 'C. elegans WS241',
    'EColiMG1655' : 'E. coli str. K-12 substr. MG1655',
    'SaureusUSA300' : 'S. aureus subsp. aureus USA300_TCH1516',
    'ScerevisiaeNCBI' : 'Saccharomyces cerevisiae S288C'
    
]
genomeTable = [ 
    'HPylori'  :  // https://www.ncbi.nlm.nih.gov/genome/169?genome_assembly_id=163534
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/093/185/GCF_000093185.1_ASM9318v1/GCF_000093185.1_ASM9318v1_genomic.fna.gz',
    'Zmobilis'     : // https://www.ncbi.nlm.nih.gov/genome/898?genome_assembly_id=300425
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/105/GCF_000007105.1_ASM710v1/GCF_000007105.1_ASM710v1_genomic.fna.gz',
    'EColiMG1655' : // https://www.ncbi.nlm.nih.gov/genome/167?genome_assembly_id=161521
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz',
    'SaureusUSA300' :  // https://www.ncbi.nlm.nih.gov/genome/154?genome_assembly_id=299284
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/017/085/GCF_000017085.1_ASM1708v1/GCF_000017085.1_ASM1708v1_genomic.fna.gz'
    'ScerevisiaeNCBI' : // https://www.ncbi.nlm.nih.gov/genome/15?genome_assembly_id=22535
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz',
    'CelegansWS241'  :
        'ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/sequence/genomic/c_elegans.PRJNA13758.WS241.genomic.fa.gz'
]
exptTable = [
    'HPylori'       : ['SRR023794,SRR023796'],
    'Zmobilis'      : ['SRR017972,SRR029606'],
    'CelegansWS241' : ['SRR065390'],
    'EColiMG1655'   : ['ERR008613'],
    'SaureusUSA300' : ['SRR022866'],
    'ScerevisiaeNCBI'   : ['SRR352384']
]

paramsKarect {
    'HPylori' : ' -matchtype=edit -celltype=haploid ',
    'Zmobilis' : ' -matchtype=edit -celltype=haploid ',
    'CelegansWS241' : ' -matchtype=hamming -celltype=diploid ',
    'EColiMG1655' : ' -matchtype=hamming -celltype=haploid ',
    'SaureusUSA300' : ' -matchtype=hamming -celltype=haploid ',
    'ScerevisiaeNCBI' : ' -matchtype=hamming -celltype=haploid '
}

paramsKarectAlign {
    'HPylori' : ' -matchtype=edit ',
    'Zmobilis' : ' -matchtype=edit ',
    'CelegansWS241' : ' -matchtype=hamming ',
    'EColiMG1655' : ' -matchtype=hamming ',
    'SaureusUSA300' : ' -matchtype=hamming ',
    'ScerevisiaeNCBI' : ' -matchtype=hamming '
}

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

exptChan = idxChan.flatMap {
    orgId, orgDesc, gnmFile -> exptTable[orgId].collect {
        [orgId, orgDesc, gnmFile, it]
    }
}.flatMap{
    orgId, orgDesc, gnmFile, exptId  ->
        sraIds = parseExptID(exptId, ',')
        sraIds.collect{
            [orgId, orgDesc, gnmFile, exptId, it]
        }
}
// .subscribe{
//     println it
// }

process SRAFetch {
    tag { orgId.toString() + " > " + exptId.toString() + " > " + sraId.toString() }

    input:
    set orgId, orgDesc, gnmFile, exptId, sraId from exptChan

    output:
    set orgId, orgDesc, gnmFile, exptId, sraId, file("${sraId}*.fastq") into pseqChan

    """
    prefetch ${sraId}
    vdb-validate '${sraId}'
    fastq-dump -I --split-files --gzip '${sraId}'
    """
}

process ConcatPariedEndFiles{
    tag { orgId.toString() + " > " + exptId.toString() + " > " + sraId.toString() }

    input:
    set orgId, orgDesc, gnmFile, exptId, sraId, file(pairedFiles) from pseqChan

    output:
    set orgId, orgDesc, gnmFile, exptId, sraId, file("${sraId}.fastq") into fseqChan
 
    """
    gunzip -c ${pairedFiles} > ${sraId}.fastq
    """
}

oexptChan = fseqChan.map{ 
    orgId, orgDesc, gnmFile, exptId, sraId, sraFile -> 
        [orgId.toString() + "-" + exptId.toString(),
         orgId, orgDesc, gnmFile, exptId, sraId, sraFile] 
}
.groupTuple()
.map{
    orgExptId, orgId, orgDesc, gnmFile, exptId, sraIds, sraFiles -> 
        [orgExptId, orgId[0], orgDesc[0], gnmFile[0],
         idxFiles[0], exptId[0], sraIds, sraFiles]
}

// oexptChan.subscribe{
//     println it
// }



process ConcatSRAFiles{
    tag { orgId.toString() + ">" + exptId.toString() }
    
    input:
    set orgExptId, orgId, orgDesc, gnmFile, exptId, sraIds, file(sraFiles) from oexptChan

    output:
    set orgExptId, orgId, orgDesc, gnmFile, exptId, sraIds, file("beforeEC.fastq") into cseqChan

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

process Karect{
    tag { orgId.toString() + ">" + exptId.toString() }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, exptId, sraIds, file(beforeEC) from beforeChan1

    output:
    set orgExptId, orgId, orgDesc, gnmFile, exptId, sraIds, file("karect_.fastq") into ecChan

    """
    karect -correct -threads=${params.nthreads} ${paramsKarect[exptId]} -inputfile=${beforeEC} 
    """ 
}


process KarectAlignBefore{
    tag {  orgId.toString() + ">" + exptId.toString() }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, exptId, sraIds, file(beforeEC) from beforeChan2

    output:
    set orgExptId, orgId, orgDesc, gnmFile, exptId, sraIds, file(beforeEC), file("beforeEC.txt") into beforeSAMChan

    """
    karect -align -threads=${params.nthreads}  ${paramsKarectAlign[exptId]}  -refgenomefile=$gnmFile  -inputfile=$beforeEC -alignfile=beforeEC.txt
    """ 
}

mergedSAMChan = beforeSAMChan
    .join(ecChan)
    .map {
        orgExptId,
            orgId1, orgDesc1, gnmFile1, idxFiles1, exptId1, sraIds1, beforeEC, beforeSAM,
            orgId2, orgDesc2, gnmFile2, idxFiles2, exptId2, sraIds2, afterEC ->
        [ orgExptId1, orgId1, orgDesc1, gnmFile1, idxFiles1,
          exptId1, sraIds1, beforeEC, beforeSAM, afterEC  ]
    }

processEvalEC{
    tag{ orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, exptId, sraIds, file(beforeEC), file(beforeSAM), file(afterEC) from mergedSAMChan

    output:
    file("result.txt") into result_channel1

// -eval -threads=${num threads} -matchtype=${match type} 
// -inputfile=${input file} -resultfile=${result file} -refgenomefile=${ref file} \
// -alignfile=${align file} -evalfile=${eval file}
    """
    karect -eval -threads=${params.nthreads}  ${paramsKarectAlign[exptId]}  -inputfile=$beforeEC  -resultfile=$afterEC -refgenomefile=$gnmFile -alignfile=$beforeSAM -evalfile=result.txt
    echo "DATASET : " $exptId  " ORGANISM : " $orgDesc  >> result.txt
    """

}


result_channel1.map{
    it.text
}.collectFile(name: 'karect_eval_sra.txt',
              storeDir: "${workflow.projectDir}",
              newLine: false)
