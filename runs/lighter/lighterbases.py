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
                aln1 = samutils.getSAMAlignment(read1.query_alignment_sequence,read1.cigarstring,read1.get_tag("MD"))
                for a,b in zip(aln1[0],aln1[1]):		#loops through and counts how many bases are aligned
			        if a == b:
                        beforeCount = beforeCount + 1.0
                
    for read2 in samfile2:
        fl2 = read2.flag
        if(fl2 & (1 << 8) != (1 << 8)):            #checks to see if this is the best version of the read
            if(fl2 & (1 << 2) != (1 << 2)):        #checks to see if the read is mapped
                aln2 = samutils.getSAMAlignment(read2.query_alignment_sequence,read2.cigarstring,read2.get_tag("MD"))
                for c,d in zip(aln2[0],aln2[1]):		#loops through and counts how many bases are aligned
			        if c == d:
                        afterCount = afterCount + 1.0
                
    print("% Increase: ", (afterCount - beforeCount)/total)
    
if __name__ == '__main__':
    stats(sys.argv[1], sys.argv[2])
                
