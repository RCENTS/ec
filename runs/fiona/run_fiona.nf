orgTable = [
    'Bpertussis18323' : 'B. pertussis 18323',
    'EcoliK12MG1655'  : 'E. coli K-12 MG1655',
    'EcoliO104H4'     : 'E. coli O104:H4 str. 2011C-3493',
    'HSapiensGRCh38'  : 'H. sapiens (GRCh38)',
    'PsyringaeB728a'  : 'P. syringae pv. syringae B728a',
    'Pfalciparum3D7'  : 'P. falciparum 3D7',
    'Scerevisae'      : 'S. cerevisae S288C',
    'SauresLGA251'    : 'S. aureus LGA251',
    'Celegans'        : 'C. elegans',
    'Dmelanogaster'   : 'D. melanogaster'
]

genomeTable = [ 
    'Bpertussis18323' : // https://www.ncbi.nlm.nih.gov/genome/1008?genome_assembly_id=170245
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/306/945/GCF_000306945.1_ASM30694v1/GCF_000306945.1_ASM30694v1_genomic.fna.gz',
    'EcoliK12MG1655' : // https://www.ncbi.nlm.nih.gov/genome/167?genome_assembly_id=161521
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz', 
    'EcoliO104H4' : // https://www.ncbi.nlm.nih.gov/genome/167?genome_assembly_id=161525
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/299/455/GCF_000299455.1_ASM29945v1/GCF_000299455.1_ASM29945v1_genomic.fna.gz',
    'HSapiensGRCh38' : // https://www.ncbi.nlm.nih.gov/genome/51?genome_assembly_id=368248
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.fna.gz',
    'PsyringaeB728a' : // https://www.ncbi.nlm.nih.gov/genome/185?genome_assembly_id=299898
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/012/245/GCF_000012245.1_ASM1224v1/GCF_000012245.1_ASM1224v1_genomic.fna.gz',
    'Pfalciparum3D7' : // https://www.ncbi.nlm.nih.gov/genome/33?genome_assembly_id=369845
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/765/GCF_000002765.4_ASM276v2/GCF_000002765.4_ASM276v2_genomic.fna.gz',
    'Scerevisae' : // https://www.ncbi.nlm.nih.gov/genome/15?genome_assembly_id=22535
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz',
    'SauresLGA251' : // https://www.ncbi.nlm.nih.gov/genome/154?genome_assembly_id=155332
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/237/265/GCF_000237265.1_ASM23726v1/GCF_000237265.1_ASM23726v1_genomic.fna.gz', 
    'Celegans'   : // 
        'ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/sequence/genomic/c_elegans.WS222.genomic.fa.gz',
    'Dmelanogaster' :
        'ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/fasta/dmel-all-chromosome-r6.24.fasta.gz'
]

exptTable = [
    // 'Bpertussis18323' : [
    //     'ERR161541'
    // ],
    'EcoliK12MG1655'  : [
        'SRR620425',
        'SRR000868',
        'ERR039477',
        'SRR611140',
        'ERR022075'
    ]
    // ,
    // 'EcoliO104H4'     : [
    //     'SRR254209'
    // ],
    // 'HSapiensGRCh38'    : [
    //     'SRR1238539'
    // ],
    // 'PsyringaeB728a'    : [
    //     'ERR005143'
    // ],
    // 'Pfalciparum3D7'  : [
    //     'ERR161543'
    // ],
    // 'Scerevisae'      : [
    //     'SRR031259',
    //     'SRR096469,SRR096470' // SRX039441
    // ],
    // 'SauresLGA251'    : [
    //     'ERR236069',
    //     'SRR070596'
    // ],
    // 'Celegans'        : [
    //     'SRR443373'
    // ],
    // 'Dmelanogaster'   : [
    //     'SRR492060',
    //     'SRR034841,SRR034842,SRR034843,SRR034844' // SRX016210
    // ]
]

runFiona = {
    'ERR161541' : 'fiona -g 4043846 '
    'ERR022075' : 'fiona_illumina -g 4639675 ',
    'SRR000868' : 'fiona -g 4639675 ',
    'ERR039477' : 'fiona -g 4639675 ',
    'SRR611140' : 'fiona -g 4639675 ',
    'SRR620425' : 'fiona -g 4639675 ',
    'SRR254209' : 'fiona -g 5273097 ',
    'SRR1238539' : 'fiona -g 2861343787 ',
    'ERR005143' : 'fiona_illumina -g 6093698 ',
    'ERR161543' : 'fiona -g  23264338 '
    'SRR031259' : 'fiona_illumina -g 12156676 '
    'SRR096469,SRR096470' : 'fiona -g 12156676 ', // SRX039441
    'ERR236069' : 'fiona -g 2799725 ',
    'SRR070596' : 'fiona -g 2799725 ',
    'SRR443373' : 'fiona_illumina -g 100286070 ',
    'SRR492060' : 'fiona_illumina -g 120381546 ',
    'SRR034841,SRR034842,SRR034843,SRR034844' : : 'fiona -g 120381546 ' // SRX016210

}

runBWA = {
    'ERR161541' : 'bwa bwasw '
    'ERR022075' : 'bwa mem ',
    'SRR000868' : 'bwa bwasw ',
    'ERR039477' : 'bwa bwasw ',
    'SRR611140' : 'bwa bwasw ',
    'SRR620425' : 'bwa bwasw ',
    'SRR254209' : 'bwa bwasw ',
    'SRR1238539' : 'bwa bwasw ',
    'ERR005143' : 'bwa mem  ',
    'ERR161543' : 'bwa bwasw '
    'SRR031259' : 'bwa mem  '
    'SRR096469,SRR096470' : 'bwa bwasw ', // SRX039441
    'ERR236069' : 'bwa bwasw ',
    'SRR070596' : 'bwa bwasw ',
    'SRR443373' : 'bwa mem  ',
    'SRR492060' : 'bwa mem  ',
    'SRR034841,SRR034842,SRR034843,SRR034844' : : 'bwa bwasw ' // SRX016210

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
    set orgId, orgDesc, gnmFile, idxFiles, exptId, sraId, file("${sraId}.fastq.gz") into fseqChan
 
    """
    gunzip -c ${pairedFiles} | gzip -1  > ${sraId}.fastq.gz
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

process ConcatSRAFiles{
    tag { orgId.toString() + " > " + exptId.toString() }
    
    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(sraFiles) from oexptChan

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file("beforeEC.fastq") into cseqChan

    script:
    """
    gunzip -c ${sraFiles} > beforeEC.fastq
    """
}

(beforeChan1, beforeChan2) = cseqChan.into(2)

process Fiona{
    tag { orgId.toString() + " > " + exptId.toString() }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC) from beforeChan1

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file("afterEC.fastq") into ecChan

    """
    ${runFiona[exptId]} -nt ${params.nthreads} ${beforeEC} afterEC.fastq
    """ 
}

process BWABefore{
    tag { orgId.toString() + " > " + exptId.toString() }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC) from beforeChan2

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file("beforeEC.sam") into beforeSAMChan

    """
    ${runBWA[exptId]} -t ${params.nthreads} ${params.genomedir}/bwa/${orgId}.fa ${beforeEC} | samtools view -bSh -F 0x900 - > bx.bam
    samtools sort -T bx.sorted -n -@ ${params.nthreads} -o beforeEC.bam bx.bam
    rm -rf bx.bam bx.sorted*
    """ 
}

process BWAAfter{
    tag { orgId.toString() + " > " + exptId.toString() }

    input:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(beforeEC), file(afterEC) from ecChan

    output:
    set orgExptId, orgId, orgDesc, gnmFile, idxFiles, exptId, sraIds, file(afterEC), file("afterEC.sam") into afterSAMChan

    """
    ${runBWA[exptId]} -t ${params.nthreads} ${params.genomedir}/bwa/${orgId}.fa ${afterEC} | samtools view -bSh -F 0x900 - > ax.bam
    samtools sort -T ax.sorted -n -@ ${params.nthreads} -o afterEC.bam ax.bam
    rm -rf ax.bam ax.sorted*
    samtools view -h -o afterEC.sam afterEC.bam
    rm -rf afterEC.bam
    """ 
}

mergedSAMChan = beforeSAMChan
    .join(afterSAMChan)
    .map {
        orgExptId,
           orgId1, orgDesc1, gnmFile1, idxFiles1, exptId1, sraIds1, beforeEC, beforeSAM,
           orgId2, orgDesc2, gnmFile2, idxFiles2, exptId2, sraIds2, afterEC, afterSAM ->
        [ orgExptId1, orgId1, orgDesc1, gnmFile1, idxFiles1,
          exptId1, sraIds1, beforeEC, afterEC, beforeSAM, afterSAM  ]
    }

processEvalEC{
    tag { orgId.toString() + " > " + exptId.toString() }

    input:
    set orgExptId, orgId, orgDesc, file(gnmFile), idxFiles, exptId, sraIds,
        file(beforeEC), file(afterEC), file(beforeSAM), file(afterSAM) from mergedSAMChan

    output:
    file("result1") into result_channel1

    """
    compute_gain -g $gnmFile --pre $beforeSAM --post $afterSAM  > result2
    """

}

result_channel1.map{
    it.text
}.collectFile(name: 'fiona_gain.txt',
              storeDir: "${workflow.projectDir}",
              newLine: false)
