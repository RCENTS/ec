import pysam
import sys

def eval_record(read_before, read_after):
    rtp = 0
    rtn = 0
    rfp = 0
    rfn = 0
    if(read_before.flag & (1 << 8) == (1 << 8)):
        return (0, 0, 0, 0, 0)
    if(read_after.flag & (1 << 8) == (1 << 8)):
        return (0, 0, 0, 0, 0)
    if(len(read_before.query_sequence) != len(read_after.query_sequence)):
        return (0, 0, 0, 0, 1)
    before_align = {}
    after_align = {}
    if(read_before.flag & (1 << 2) != (1 << 2)):
        for x,y,z in read_before.get_aligned_pairs(with_seq=True):
            if x:
                before_align[x] = (y, z)
        # samutils.getSAMAlignment(read1.query_alignment_sequence,read1.cigarstring,read1.get_tag("MD"))
    if(read_after.flag & (1 << 2) != (1 << 2)):
        for x,y,z in read_after.get_aligned_pairs(with_seq=True):
            if x:
                after_align[x] = (y, z)
        #samutils.getSAMAlignment(read2.query_alignment_sequence,read2.cigarstring,read2.get_tag("MD"))
    for pos in set(before_align.keys() + after_align.keys()):
        beforeMapped = False
        afterMapped = False
        if pos in before_align and before_align[pos][0]:
            beforeMapped = True
        if pos in after_align and after_align[pos][0]:
            afterMapped = True
        if (beforeMapped and not afterMapped):
            rfp = rfp + 1
        if (not beforeMapped and afterMapped):
            rtp = rtp + 1
        if (not beforeMapped and not afterMapped):
            rfn = rfn + 1
        if (beforeMapped and afterMapped):
            rtn = rtn + 1     
    return(rtp, rtn, rfp, rfn, 0)


def stats(before_sam, after_sam):
    samfile1 = pysam.AlignmentFile(before_sam, "rb")
    samfile2 = pysam.AlignmentFile(after_sam, "rb")    #converting bam to sam
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    df = 0
    rid1 = ''
    rid2 = ''
    try:
        while True:
            rx1 = samfile1.next()
            rid1 = rx1.query_name
            while rid2 != rid1:
                rx2 = samfile2.next()
                rid2 = rx2.query_name
            if rid1 == rid2:
                rtp, rfp, rtn, rfn, rdx = eval_record(rx1, rx2)
                tp = tp + rtp
                fp = fp + rfp
                tn = tn + rtn
                fn = fn + rfn
                df = df + rdx 
    except StopIteration as e:
        pass
    except IOError as e:
        print str(e)
    samfile1.close()
    samfile2.close()
    #
    # print("tp: ", tp)
    # print("fp: ", fp)
    # print("tn: ", tn)
    # print("fn: ", fn)
    #print("count: ", count)
    sensitivity = float(tp)/float((tp+fn))
    specificity = float(tn)/float((tn+fp))
    precision = float(tp)/float((tp+fp))
    gain = float((tp-fp))/float((tp+fn))
    #
    # print("sensitivity: ", sensitivity)
    # print("specificity: ", specificity)
    # print("precision: ", precision)
    # print("gains: ", gain)
    return [tp, fp, tn, fn, df,
            sensitivity, specificity, precision, gain]

def main(before_sam, after_sam, run_id):
    eval_result = stats(before_sam, after_sam)
    print "\t".join([run_id] + [str(x) for x in eval_result])

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
