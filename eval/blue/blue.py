import pysam
import samutils
import sys

def stats(sam):
    x = same
    samfile = pysam.AlignmentFile(x, "rb")
    edits = [0,0,0,0,0,0,0,0,0,0,0,0]
    c = 0
    for read in samfile:
        		#print "read id : ", read.query_name
		    #print "read    : ", read.query_alignment_sequence
		    #print "cigar   : ", read.cigarstring
		    #print "tags    : ", read.get_tag("MD")
		    #print read
		    fl = read.flag
		    if(fl & (1 << 8) != (1 << 8)):			#checks to see if this is the best version of the read
			    if(fl & (1 << 2) == (1 << 2)):		#checks to see if the read is unmapped
				    edits[11] = edits[11] + 1	
			    else:
				    aln = samutils.getSAMAlignment(read.query_alignment_sequence,read.cigarstring,read.get_tag("MD"))
				    #print aln[0]
				    #print aln[1]
				    count = 0
				    for a,b in zip(aln[0],aln[1]):		#loops through and counts how many edits need to be made
					    if a != b:
						    count = count + 1
				    if(count >= 10):
					    edits[10] = edits[10] + 1
				    else:
					    edits[count] = edits[count] + 1
    print edits
    return edits
    
if __name__ == '__main__':
    stats(sys.argv[1])
