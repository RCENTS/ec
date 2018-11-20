orgTable = [
    'EColiMG1655' : 'E. coli K-12 MG1655',
    'MtuberculosisH37Rv' : 'M. tuberculosis H37Rv',
    'SentericaSL476' : 'S. enterica subsp. enterica serovar Heidelberg str. SL476',
    'LmonocytogenesFSLR2561' : 'L. monocytogenes FSL R2-561',
    'PsyringaeB728a' : 'P. syringae pv. syringae B728a',
    'BdentiumBd1' : 'B. dentium Bd1',
    'OtsutsugamushiBor' : 'O. tsutsugamushi str. Boryong'
]

genomeTable = [
    'EColiMG1655' : 
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz',
    'MtuberculosisH37Rv' :  // https://www.ncbi.nlm.nih.gov/genome/166?genome_assembly_id=159857
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz',
    'SentericaSL476' : // https://www.ncbi.nlm.nih.gov/genome/152?genome_assembly_id=299234
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/020/705/GCF_000020705.1_ASM2070v1/GCF_000020705.1_ASM2070v1_genomic.fna.gz',
    'LmonocytogenesFSLR2561' : // https://www.ncbi.nlm.nih.gov/genome/159?genome_assembly_id=159666
      'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/168/575/GCF_000168575.2_ASM16857v2/GCF_000168575.2_ASM16857v2_genomic.fna.gz',
    'PsyringaeB728a' : // https://www.ncbi.nlm.nih.gov/genome/185?genome_assembly_id=299898
      'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/012/245/GCF_000012245.1_ASM1224v1/GCF_000012245.1_ASM1224v1_genomic.fna.gz',
    'BdentiumBd1' : // https://www.ncbi.nlm.nih.gov/genome/809?genome_assembly_id=168861
      'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/024/445/GCF_000024445.1_ASM2444v1/GCF_000024445.1_ASM2444v1_genomic.fna.gz',
    'OtsutsugamushiBor' : // https://www.ncbi.nlm.nih.gov/genome/710?genome_assembly_id=168255
      'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/063/545/GCF_000063545.1_ASM6354v1/GCF_000063545.1_ASM6354v1_genomic.fna.gz'
]

genomeNumTable = [
    'EColiMG1655' : '4641650',
    'MtuberculosisH37Rv' : '4411532',
    'SentericaSL476' : '4888768',
    'LmonocytogenesFSLR2561' : '2973801',
    'PsyringaeB728a' : '6093698',
    'BdentiumBd1' : '2636367',
    'OtsutsugamushiBor' : '2127051'
]

exptTable = [
    'EColiMG1655' : ['SRR519926', 
                     'SRR522163'],
    'MtuberculosisH37Rv' : ['SRR1200797'],
    'SentericaSL476' : ['SRR1203044',
                        'SRR1206093'],
    'LmonocytogenesFSLR2561' : ['SRR1198952'],
    'PsyringaeB728a' : ['SRR1119292'],
    'BdentiumBd1' : ['SRR1151311'],
    'OtsutsugamushiBor' : ['SRR1202083']
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

exptChan = idxChan.flatMap {
    orgId, orgDesc, gnmFile -> exptTable[orgId].collect {
        [orgId, orgDesc, gnmFile, it]
    }
}.flatMap{
    orgId, orgDesc, gnmFile, exptId  ->
        sraIds = parseExptID(exptId, ',')
        sraIds.collect{
            [orgId, orgDesc, gnmFile,exptId, it]
        }
}
//.subscribe{
//     println it
//}

process SRAFetch {
    tag { orgDesc.toString() + " > " + exptId.toString() + " > " + sraId.toString() }

    input:
    set orgId, orgDesc, gnmFile, exptId, sraId from exptChan

    output:
    set orgId, orgDesc, gnmFile, exptId, sraId, file("${sraId}*.fastq.gz") into pseqChan

    """
    prefetch ${sraId}
    vdb-validate '${sraId}'
    fastq-dump -I --split-files --gzip '${sraId}'
    """
}

process ConcatPariedEndFiles{
    tag { orgDesc.toString() + " > " + exptId.toString() + " > " + sraId.toString() }

    input:
    set orgId, orgDesc, gnmFile, exptId, sraId, file(pairedFiles) from pseqChan

    output:
    set orgId, orgDesc, gnmFile, exptId, sraId, file("${sraId}.fastq.gz") into fseqChan
 
    """
    gunzip -c ${pairedFiles}  | gzip -1 > ${sraId}.fastq.gz
    """
}

oexptChan = fseqChan.map{ 
    orgId, orgDesc, gnmFile, exptId, sraId, sraFile -> 
        [orgDesc.toString() + "-" + exptId.toString(),
         orgId, orgDesc, gnmFile, exptId, sraId, sraFile] 
}
.groupTuple()
.map{
    orgExptId, orgId, orgDesc, gnmFile, exptId, sraIds, sraFiles -> 
        [orgExptId, orgId[0], orgDesc[0], gnmFile[0],
         exptId[0], sraIds, sraFiles]
}

//oexptChan.subscribe{
//     println it
//}

process ConcatSRAFiles{
    tag { orgDesc.toString() + " > " + exptId.toString() }
    
    input:
    set orgExptId, orgId, orgDesc, gnmFile, exptId, sraIds, file(sraFiles) from oexptChan

    output:
    set orgExptId, orgId, orgDesc, gnmFile, exptId, sraIds, file("beforeEC.fastq") into cseqChan

    """
    gunzip -c ${sraFiles} > beforeEC.fastq
    """
}

process ACE{
    tag { orgDesc.toString() + " > " + exptId.toString() }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, exptId, sraIds, file(beforeEC) from cseqChan

    output:
    set orgExptId, orgId, orgDesc, gnmFile, exptId, sraIds, file(beforeEC), file("afterEC.fastq") into ecChan

    """
    ACE ${genomeNumTable[orgId]} ${beforeEC} afterEC.fastq ${params.lognthreads}
    """ 
}

(evalChan1, evalChan2) = ecChan.into(2)


process EvalECReads {
    tag{ orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, exptId, sraIds, file(beforeEC), file(afterEC) from evalChan1

    output:
    file(result1) into result_channel1

    """
    echo "DATASET : " $exptId  " ORGANISM : " $orgDesc  > result1
    readSearch $gnmFile $beforeEC $afterEC >> result1
    echo "-----END-----" >> result1
    """ 

}


process EvalECReads {
    tag{ orgExptId.replace('-SRR', ' > SRR') }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, exptId, sraIds, file(beforeEC), file(afterEC) from evalChan2

    output:
    file(result2) into result_channel2

    """
    echo "DATASET : " $exptId  " ORGANISM : " $orgDesc  > result2
    kmerSearch ${params.kmer} $gnmFile $beforeEC $afterEC >> result2
    echo "-----END-----" >> result2
    """ 

}


result_channel1.map{
    it.text
}.collectFile(name: 'ace_eval_reads_miseq.txt',
              storeDir: "${workflow.projectDir}",
              newLine: false)

result_channel2.map{
    it.text
}.collectFile(name: 'ace_eval_bases_miseq.txt',
              storeDir: "${workflow.projectDir}",
              newLine: false)

