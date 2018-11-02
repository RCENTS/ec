import pysam
import sys


def zdiv(x, nr):
    if nr > 0:
        return float(x)/float(nr)
    else:
        return 0.0

def eval_record(read_before, read_after):
    before = True
    after = True
    fl1 = read_before.flag
    fl2 = read_after.flag
    if(fl1 & (1 << 8) != (1 << 8)) and (fl2 & (1 << 8) != (1 << 8)):  #primary alignment
        if(fl1 & (1 << 2) == (1 << 2)):        #checks to see if the read is unmapped
            before = False
        if(fl2 & (1 << 2) == (1 << 2)):
            after = False
        if (not before) and after:
            return (1, 0, 0, 0) # TP
        if before and (not after):
            return (0, 1, 0, 0) # FP
        if before and after:
            return (0, 0, 1, 0) # TN
        if (not before) and (not after):
            return (0, 0, 0, 1) # FN
    return (0, 0, 0, 0)


def stats(before_sam, after_sam):
    samfile1 = pysam.AlignmentFile(before_sam, "rb")
    samfile2 = pysam.AlignmentFile(after_sam, "rb")    #converting bam to sam
    tp = 0
    fp = 0
    tn = 0
    fn = 0
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
                rtp, rfp, rtn, rfn = eval_record(rx1, rx2)
                tp = tp + rtp
                fp = fp + rfp
                tn = tn + rtn
                fn = fn + rfn
    except StopIteration as e:
        pass
    except IOError as e:
        print str(e)
    samfile1.close()
    samfile2.close()
    #
    print("tp: ", tp)
    print("fp: ", fp)
    print("tn: ", tn)
    print("fn: ", fn)
    #print("count: ", count)
    sensitivity = zdiv(float(tp), float((tp+fn)))
    specificity = zdiv(float(tn), float((tn+fp)))
    precision = zdiv(float(tp), float((tp+fp)))
    gain = zdiv(float((tp-fp)), float((tp+fn)))
    #
    # print("sensitivity: ", sensitivity)
    # print("specificity: ", specificity)
    # print("precision: ", precision)
    # print("gains: ", gain)
    return [tp, fp, tn, fn, 
            sensitivity, specificity, precision, gain]

def main(before_sam, after_sam, run_id):
    eval_result = stats(before_sam, after_sam)
    print "\t".join([run_id] + [str(x) for x in eval_result])

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
