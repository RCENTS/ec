
orgTable = [
    'SaureusMW2'  : 'Staphylococcus aureus',
    'Hacinonychis':'Helicobacter acinonychis str. Sheeba'
]

genomeTable = [ 
    'SaureusMW2' :
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/265/GCF_000011265.1_ASM1126v1/GCF_000011265.1_ASM1126v1_genomic.fna.gz', 
    'Hacinonychis': //https://www.ncbi.nlm.nih.gov/genome/1444?genome_assembly_id=300646
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/305/GCF_000009305.1_ASM930v1/GCF_000009305.1_ASM930v1_genomic.fna.gz'
]

exptTable = [
    // 'Hacinonychis':'http://sharcgs.molgen.mpg.de/data/reads_seq.gz',
    'SaureusMW2' : 'http://www.genomic.ch/edena/mw2Reads.seq.gz'
]


paramsTable = [
    'SaureusMW2' : '  2820460  1.00 ',
    'Hacinonychis' : '  1557590   1.6 '
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


process SeqDownload {
    tag{ orgDesc }

    input:
    set orgId, orgDesc, gnmFile from orgChan

    output:
    set orgId, orgDesc, gnmFile, file("${orgId}.fasta") into pseqChan

    """
    wget ${exptTable[orgId]} -O ${orgId}.fasta.gz
    gunzip ${orgId}.fasta.gz
    """
}

process HiTEC{
    tag{ orgDesc }

    input:
    set orgId, orgDesc, gnmFile, file(beforeEC) from pseqChan

    output:
    set orgId, orgDesc, gnmFile, file(beforeEC), file("afterEC.fasta") into ecChan

    """
    hitec ${beforeEC} afterEC.fasta ${paramsTable[orgId]}
    """ 
}

process EvalEC{
    tag{ orgDesc }

    input:
    set orgId, orgDesc, gnmFile, file(beforeEC), file(afterEC) from ecChan

    output:
    file(result) into result_channel

    """
    echo  "DATASET : " $orgDesc > result
    readSearch $gnmFile $beforeEC $afterEC  >> result
    """ 

}


result_channel.map{
    it.text
}.collectFile(name: 'hitec_online_gain.txt',
              storeDir: "${workflow.projectDir}",
              newLine: false)
