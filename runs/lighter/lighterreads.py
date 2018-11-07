import pysam
import sys

def zdiv(x, nr):
    if nr > 0:
        return float(x)/float(nr)
    else:
        return 0.0

def stats(before_sam, after_sam):
    samfile1 = pysam.AlignmentFile(before_sam, "rb")
    samfile2 = pysam.AlignmentFile(after_sam, "rb")    #converting bam to sam

    beforeCount = 0
    afterCount = 0
    totalReads = 0
    for read1 in samfile1:
        fl1 = read1.flag
        if(fl1 & (1 << 8) != (1 << 8)):            #checks to see if this is the best version of the read
            totalReads += 1
            if(fl1 & (1 << 2) != (1 << 2)):        #checks to see if the read is mapped
                beforeCount = beforeCount + 1.0
                
    for read2 in samfile2:
        fl2 = read2.flag
        if(fl2 & (1 << 8) != (1 << 8)):            #checks to see if this is the best version of the read
            if(fl2 & (1 << 2) != (1 << 2)):        #checks to see if the read is mapped
                afterCount = afterCount + 1.0

    samfile1.close()
    samfile2.close()
    return [beforeCount, afterCount, totalReads,
            zdiv((afterCount - beforeCount), totalReads)]
    #print("% Increase: ", (afterCount - beforeCount)/total)
    
if __name__ == '__main__':
    lst = stats(sys.argv[0], sys.argv[1])
    print "\t".join([sys.argv[2]] + [str(x) for x in lst])
                
    
