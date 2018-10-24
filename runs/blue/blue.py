import pysam
import samutils as su
import sys

def zdiv(x, nr):
    if nr > 0:
        return float(x)/float(nr)
    else:
        return 0.0

def sam_stats(inFile):
    samfile = pysam.AlignmentFile(inFile, "rb")
    edits = [0,0,0,0,0,0,0,0,0,0,0,0]
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
                edits[11] = edits[11] + 1    
            else:
                aln = su.getSAMAlignment(read.query_alignment_sequence,read.cigarstring,read.get_tag("MD"))
                #print aln[0]
                #print aln[1]
                count = 0
                for a,b in zip(aln[0],aln[1]):        #loops through and counts how many edits need to be made
                    if a != b:
                        count = count + 1
                if(count >= 10):
                    edits[10] = edits[10] + 1
                else:
                    edits[count] = edits[count] + 1
    #print edits
    return [nreads] + edits + [zdiv(x, nreads) for x in edits]

def main(beforeEC, afterEC, oID):
    bstats = sam_stats(beforeEC)
    astats = sam_stats(afterEC)
    bstats_str = "\t".join([oID, "BLUE", "BEFORE"] + [str(x) for x in bstats])
    astats_str = "\t".join([oID, "BLUE", "AFTER"] + [str(x) for x in astats])
    print bstats_str
    print astats_str

    
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
