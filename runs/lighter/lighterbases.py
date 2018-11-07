import pysam
import sys

def zdiv(x, nr):
    if nr > 0:
        return float(x)/float(nr)
    else:
        return 0.0

def eval_record(read):
    aln_bases = 0
    err_bases = 0
    if(read.flag & (1 << 8) != (1 << 8)):
        if(read.flag & (1 << 2) == (1 << 2)): # read is unmapped!
            aln_bases = 0
            err_bases = len(read.query_sequence) #all bases in the read
        else:
            #print read1.query_sequence
            rbp = read.get_aligned_pairs(matches_only=True, with_seq=True)
            for _, _, c in rbp:
                if c.isupper():
                    aln_bases = aln_bases + 1  #Number of perfect bases the before sequence
                else:
                    err_bases = err_bases + 1 # no. of mismatch bases
            cgx = read.get_cigar_stats()[0]
            err_bases = err_bases + cgx[1] # insertions
            err_bases = err_bases + cgx[2] # deletions
    return (aln_bases, err_bases)


def stats(before_sam, after_sam):
    samfile1 = pysam.AlignmentFile(before_sam, "rb")
    samfile2 = pysam.AlignmentFile(after_sam, "rb")    #converting bam to sam
    beforeCount = 0
    afterCount = 0
    totalBases = 0
    for rx1 in samfile1:
        totalBases += len(rx1.query_sequence)
        albases, _ = eval_record(rx1)
        beforeCount += albases
    for rx2 in samfile1:
        albases, _ = eval_record(rx2)
        afterCount += albases
    samfile1.close()
    samfile2.close()
    return [beforeCount, afterCount, totalBases,
            zdiv((afterCount - beforeCount), totalBases)]
    #print("% Increase: ", (afterCount - beforeCount)/total)
    
if __name__ == '__main__':
    lst = stats(sys.argv[1], sys.argv[2])
    print "\t".join([sys.argv[2]] + [str(x) for x in lst])
