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
    count = 0
    for read1 in samfile1:
        before = True
        after = True
        fl1 = read1.flag
        if(fl1 & (1 << 8) != (1 << 8)):            #checks to see if this is the best version of the read
            if(fl1 & (1 << 2) == (1 << 2)):        #checks to see if the read is unmapped
                before = False
            read1_id = read1.query_name
            try:
                read2_found = samfile2_index.find(read1_id)
            except KeyError:
                count = count + 1
                continue
            for read2 in read2_found:        #loops through all versions with the same read id
                fl2 = read2.flag
                if(fl2 & (1 << 8) != (1 << 8)):
                    if(fl2 & (1 << 2) == (1 << 2)):
                        after = False
                    if(before and not after):    #checks each condition and increments accordingly
                        fp = fp + 1.0
                    if(not before and after):
                        tp = tp + 1.0
                    if(not before and not after):
                        fn = fn + 1.0
                    if(before and after):
                        tn = tn + 1.0
                    break

    print("tp: ", tp)
    print("fp: ", fp)
    print("tn: ", tn)
    print("fn: ", fn)
    #print("count: ", count)

    print("sensitivity: ", tp/(tp+fn))
    print("specificity: ", tn/(tn+fp))
    print("precision: ", tp/(tp+fp))
    print("gains: ", (tp-fp)/(tp+fn))
    return {'TP': tp}

if __name__ == '__main__':
    stats(sys.argv[1], sys.argv[2])
