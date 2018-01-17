import pysam
import samutils
import sys

def stats(before_sam, after_sam):
    samfile1 = pysam.AlignmentFile(before_sam, "rb")
    samfile2 = pysam.AlignmentFile(after_sam, "rb")    #converting bam to sam

    beforeCount = 0.0
    afterCount = 0.0
    total = 0.0
    for read1 in samfile1:
        fl1 = read1.flag
        if(fl1 & (1 << 8) != (1 << 8)):            #checks to see if this is the best version of the read
            total = total + 1.0
            if(fl1 & (1 << 2) != (1 << 2)):        #checks to see if the read is mapped
                beforeCount = beforeCount + 1.0
                
    for read2 in samfile2:
        fl2 = read2.flag
        if(fl2 & (1 << 8) != (1 << 8)):            #checks to see if this is the best version of the read
            if(fl2 & (1 << 2) != (1 << 2)):        #checks to see if the read is mapped
                afterCount = afterCount + 1.0
                
    print("% Increase: ", (afterCount - beforeCount)/total)
    
if __name__ == '__main__':
    stats(sys.argv[1], sys.argv[2])
                
    
