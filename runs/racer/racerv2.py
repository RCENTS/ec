import pysam
import sys

def stats(before_sam, after_sam, name):
    def eval_record(read_before, read_after):
        err_before = 0.0
        err_after = 0.0
        fl1 = read_before.flag
        if(fl1 & (1 << 8) != (1 << 8)):
            if(fl1 & (1 << 2) == (1 << 2)):
                err_before = len(read_before.query_sequence) # unmapped reads = all bases are error
            else:
                rbp = read_before.get_aligned_pairs(matches_only=True, with_seq=True)
                for _,_,c in rbp:
                    if c.islower():
                        err_before = err_before + 1 # no. of mismatch bases
                cgx = read_before.get_cigar_stats()[0]
                err_before = err_before + cgx[1] # no. of insertions
                err_before = err_before + cgx[2] # no. of deletions
        else:
            return (err_after, err_before)
        fl2 = read_after.flag
        if(fl2 & (1 << 8) != (1 << 8)):
            if(fl2 & (1 << 2) == (1 << 2)):
                err_after = len(read_after.query_sequence)
            else:
                rbp = read_after.get_aligned_pairs(matches_only=True, with_seq=True)
                for _,_,c in rbp:
                    if c.islower():
                        err_after = err_after + 1 # no. of mismatch bases
                cgx = read_after.get_cigar_stats()[0]
                err_after = err_after + cgx[1] # no. of insertions
                err_after = err_after + cgx[2] # no. of deletions
        return (err_after, err_before)
    samfile1 = pysam.AlignmentFile(before_sam, "rb")
    samfile2 = pysam.AlignmentFile(after_sam, "rb")
    dset_err_before = 0.0
    dset_err_after = 0.0
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
                xea, xeb = eval_record(rx1, rx2)
                dset_err_after = dset_err_after + xea
                dset_err_before = dset_err_before + xeb
    except StopIteration as e:
        pass
    except IOError as e:
        print str(e)
    samfile1.close()
    samfile2.close()
    return ((dset_err_before-dset_err_after)/dset_err_before) if dset_err_before > 0 else 0

if __name__ == '__main__':
    org = sys.argv[3]
    print org, "RACER", stats(sys.argv[1], sys.argv[2], org)
