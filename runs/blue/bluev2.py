import pysam
import sys

def zdiv(x, nr):
    if nr > 0:
        return float(x)/float(nr)
    else:
        return 0.0

def sam_stats(inFile):
    samfile = pysam.AlignmentFile(inFile, "rb")
    edit_histo = [0,0,0,0,0,0,0,0,0,0,0,0]
    nreads = 0
    for read in samfile:
        #print "read id : ", read.query_name
        #print "read    : ", read.query_alignment_sequence
        #print "cigar   : ", read.cigarstring
        #print "tags    : ", read.get_tag("MD")
        #print read
        fl = read.flag
        if(fl & (1 << 8) != (1 << 8)):       #checks to see if this is the best version of the read
            nreads += 1
            if(fl & (1 << 2) == (1 << 2)):   #checks to see if the read is unmapped
                edit_histo[11] = edit_histo[11] + 1    
            else:
                read_error = 0
                rbp = read.get_aligned_pairs(matches_only=True, with_seq=True)
                for _,_,c in rbp:
                    if c.islower():
                        read_error = read_error + 1 # no. of mismatch bases
                cgx = read.get_cigar_stats()[0]
                read_error = read_error + cgx[1] # no. of insertions
                read_error = read_error + cgx[2] # no. of deletions
                if(read_error >= 10):
                    edit_histo[10] = edit_histo[10] + 1
                else:
                    edit_histo[read_error] = edit_histo[read_error] + 1
    #print edits
    samfile.close()
    return [nreads] + edit_histo + [zdiv(x, nreads) for x in edit_histo]

def main(beforeEC, afterEC, oID):
    bstats = sam_stats(beforeEC)
    astats = sam_stats(afterEC)
    bstats_str = "\t".join([oID, "BLUE", "BEFORE"] + [str(x) for x in bstats])
    astats_str = "\t".join([oID, "BLUE", "AFTER"] + [str(x) for x in astats])
    print bstats_str
    print astats_str

    
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
