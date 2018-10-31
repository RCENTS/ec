import pysam
import samutils as su
import sys

def stats(before_sam, after_sam, name):
    def eval_record(read1, read2):
        reb = 0.0
        rea = 0.0
        fl1 = read1.flag
        if(fl1 & (1 << 8) != (1 << 8)):
            if(fl1 & (1 << 2) == (1 << 2)):
                reb = len(read1.query_sequence) #In case of unmapped, adds in all bases in the read
            else:
                #print read1.query_sequence
                aln1 = su.getSAMAlignment(read1.query_sequence,read1.cigarstring,read1.get_tag("MD"))
                for a,b in zip(aln1[0],aln1[1]):
                    if a != b:
                        reb = reb + 1        #Number of errors in the before sequence
            # print("eb", eb)
            # print("ea", ea)
        else:
            return (rea, reb)
        fl2 = read2.flag
        if(fl2 & (1 << 8) != (1 << 8)):
            if(fl2 & (1 << 2) == (1 << 2)):
                rea = len(read2.query_sequence)
            else:
                aln2 = su.getSAMAlignment(read2.query_sequence,read2.cigarstring,read2.get_tag("MD"))
                for a,b in zip(aln2[0],aln2[1]):
                    if a != b:
                        rea = rea + 1    #Number of errors in the after sequence
        return (rea, reb)
    samfile1 = pysam.AlignmentFile(before_sam, "rb")
    samfile2 = pysam.AlignmentFile(after_sam, "rb")
    eb = 0.0
    ea = 0.0
    rid1 = ''
    rid2 = ''
    try:
        while True:
            read1 = samfile1.next()
            rid1 = read1.query_name
            while rid2 != rid1:
                read2 = samfile2.next()
                rid2 = read2.query_name
            if rid1 == rid2:
                xea, xeb = eval_record(read1, read2)
                ea = ea + xea
                eb = eb + xeb
    except StopIteration as e:
        pass
    except IOError as e:
        print str(e)
    return ((eb-ea)/eb) if eb > 0 else 0

if __name__ == '__main__':
    org = sys.argv[3]
    print org, "RACER", stats(sys.argv[1], sys.argv[2], org)
