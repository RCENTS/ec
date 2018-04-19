import pysam
import samutils
import sys

def stats(before_sam, after_sam, name):
    samfile1 = pysam.AlignmentFile(before_sam, "rb")
    samfile2 = pysam.AlignmentFile(after_sam, "rb")
    samfile2_index = pysam.IndexedReads(samfile2)
    samfile2_index.build()

    eb = 0.0
    ea = 0.0
    print("start")
    for read1 in samfile1:
        fl1 = read1.flag
        if(fl1 & (1 << 8) != (1 << 8)):
            if(fl1 & (1 << 2) == (1 << 2)):
                eb = eb + len(read1.query_sequence)    #In case of unmapped, adds in all bases in the read
            else:
                            #print read1.query_sequence
                aln1 = samutils.getSAMAlignment(read1.query_sequence,read1.cigarstring,read1.get_tag("MD"))
                for a,b in zip(aln1[0],aln1[1]):
                    if a != b: 
                        eb = eb + 1        #Number of errors in the before sequence
            print("eb", eb)
            print("ea", ea)
        else:
            continue
        read1_id = read1.query_name
        try:
            read2_found = samfile2_index.find(read1_id)
        except KeyError:
            continue
        for read2 in read2_found:
            fl2 = read2.flag
            if(fl2 & (1 << 8) != (1 << 8)):
                if(fl2 & (1 << 2) == (1 << 2)):
                    ea = ea + len(read2.query_sequence)
                else:
                    aln2 = samutils.getSAMAlignment(read2.query_sequence,read2.cigarstring,read2.get_tag("MD"))
                    for a,b in zip(aln2[0],aln2[1]):
                        if a != b:
                            ea = ea + 1    #Number of errors in the after sequence
                break

    
    # print((eb-ea)/eb)
    return((eb-ea)/eb)
    
if __name__ == '__main__':
    org = sys.argv[3]
    print org, stats(sys.argv[1], sys.argv[2], org)
