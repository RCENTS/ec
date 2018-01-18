import pysam
import samutils
import sys


def stats(before_sam, after_sam):
    samfile1 = pysam.AlignmentFile(before_sam, "rb")
    samfile2 = pysam.AlignmentFile(after_sam, "rb")    #converting bam to sam
    samfile2_index = pysam.IndexedReads(samfile2)
    samfile2_index.build()

    fp = 0.0
    tp = 0.0
    fn = 0.0
    tn = 0.0
    diflen = 0
    
    for read1 in samfile1:
        fl1 = read1.flag
        if(fl1 & (1 << 8) != (1 << 8)):
            len1 = len(read1.query_sequence)
            read1_id = read1.query_name
            try:
                read2_found = samfile2_index.find(read1_id)
            except KeyError:
                continue
            for read2 in read2_found:        #loops through all versions with the same read id
                fl2 = read2.flag
                if(fl2 & (1 << 8) != (1 << 8)):
                    len2 = len(read2.query_sequence)
                    if(len1 != len2):
                        diflen = diflen + 1
                    else:
                        read1mapped = True;
                        read2mapped = True;
                        if(fl1 & (1 << 2) == (1 << 2)):
                            read1mapped = False
                        else:
                            aln1 = samutils.getSAMAlignment(read1.query_alignment_sequence,read1.cigarstring,read1.get_tag("MD"))
                        if(fl2 & (1 << 2) == (1 << 2)):
                            read2mapped = False
                        else:
                            aln2 = samutils.getSAMAlignment(read2.query_alignment_sequence,read2.cigarstring,read2.get_tag("MD"))
                        x = 0
                        y = 0
                        while x < len1 and y < len2:
                            beforeMapped = True
                            afterMapped = True
                            while aln1[1][x] == '_':
                                x = x + 1
                            while aln2[1][y] == '_':
                                y = y + 1
                            if aln1[0][x] != aln1[1][x] or not read1mapped:
                                beforeMapped = False
                            if aln2[0][x] != aln2[1][x] or not read2mapped:
                                afterMapped = False    
                            if(beforeMapped and not afterMapped):    #checks each condition and increments accordingly
                                fp = fp + 1.0
                            if(not beforeMapped and afterMapped):
                                tp = tp + 1.0
                            if(not beforeMapped and not afterMapped):
                                fn = fn + 1.0
                            if(beforeMapped and afterMapped):
                                tn = tn + 1.0      
                            x = x + 1
                            y = y + 1
                            
    print("tp: ", tp)
    print("fp: ", fp)
    print("tn: ", tn)
    print("fn: ", fn)
   
    print("sensitivity: ", tp/(tp+fn))
    print("specificity: ", tn/(tn+fp))
    print("precision: ", tp/(tp+fp))
    print("gains: ", (tp-fp)/(tp+fn))
    return {'TP': tp}

if __name__ == '__main__':
    stats(sys.argv[1], sys.argv[2])
                            
