import pysam
import sys

def stats(before_sam, after_sam):
    def eval_record(read1, read2):
        reb = 0
        rea = 0
        fl1 = read1.flag
        bpfx = 0
        apfx = 0
        berr = 0
        aerr = 0
        if(fl1 & (1 << 8) != (1 << 8)):
            if(fl1 & (1 << 2) == (1 << 2)): # read is unmapped!
                reb = 0
                berr = len(read1.query_sequence) #all bases in the read
            else:
                #print read1.query_sequence
                rbp = read1.get_aligned_pairs(matches_only=True, with_seq=True)
                for a,b,c in rbp:
                    if c.isupper():
                        reb = reb + 1        #Number of perfect bases the before sequence
                    else:
                        berr = berr + 1 # no. of mismatch bases
                cgx = read1.get_cigar_stats()[0]
                berr = berr + cgx[1] # insertions
                berr = berr + cgx[2] # deletions
                bpfx = (1 if berr == 0 else 0) # perfect 
            # print("eb", eb)
            # print("ea", ea)
        else:
            return (0, 0, 0, 0)
        fl2 = read2.flag
        if(fl2 & (1 << 8) != (1 << 8)):
            if(fl2 & (1 << 2) == (1 << 2)): # unmapped
                rea = 0 # read2.query_sequence
            else:
                rap = read2.get_aligned_pairs(matches_only=True, with_seq=True)
                for a,b,c in rap:
                    if c.isupper():
                        rea = rea + 1    #Number of perfect bases in the after sequence
                    else:
                        aerr = aerr + 1 # no. of mismatch bases
                cgx = read2.get_cigar_stats()[0]
                aerr = aerr + cgx[1] # insertions
                aerr = aerr + cgx[2] # deletions
                apfx = (1 if aerr == 0 else 0) # perfect 
        else:
            return (0, 0, bpfx, 0)
        if rea > reb: 
            return (1, 0, bpfx, apfx)
        if rea < reb: 
            return (0, 1, bpfx, apfx)
        return (0, 0, bpfx, apfx)
    samfile1 = pysam.AlignmentFile(before_sam, "rb")
    samfile2 = pysam.AlignmentFile(after_sam, "rb")
    rbetter = 0
    rworse = 0
    rbpfx = 0
    rapfx = 0
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
                xb, xw, px0, px1 = eval_record(read1, read2)
                rbetter = rbetter + xb
                rworse = rworse + xw
                rbpfx = rbpfx + px0 
                rapfx = rapfx + px1
    except StopIteration as e:
        pass
    except IOError as e:
        print str(e)
    return [str(rbetter), str(rworse), str(rbpfx), str(rapfx)]

if __name__ == '__main__':
    org = sys.argv[3]
    print "\t".join([org, "RACER"] + stats(sys.argv[1], sys.argv[2]))
