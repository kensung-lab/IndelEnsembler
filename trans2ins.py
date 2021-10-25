from __future__ import print_function
import sys, pyfaidx

results = sys.stdin.readlines()
fasta = pyfaidx.Fasta(sys.argv[1])

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
def reverse_complement(seq):
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

class tra_bkp_t:
    def __init__(self, line):
        sl = line.split()
        self.id = sl[0]
        self.dest_chr, self.src_chr = sl[1], sl[4]
        self.dst_pos, self.src_pos = int(sl[2]), int(sl[5])
        self.dst_str, self.src_str = sl[3], sl[6]
        self.precision = sl[8]

    def is_rc(self):
        return self.src_str == self.dst_str


fwd_calls, rev_calls = dict(), dict()
for line in results:
    bp = tra_bkp_t(line)
    if bp.dst_str == '+':
        fwd_calls[bp.id] = bp
    else:
        rev_calls[bp.id] = bp

for id in fwd_calls:
    if id not in rev_calls: continue

    fwd_call, rev_call = fwd_calls[id], rev_calls[id]

    dst_chr, dst_fwd_pos, dst_rev_pos = fwd_call.dest_chr, fwd_call.dst_pos, rev_call.dst_pos
    mh_seq, mh_len = "", 0
    if dst_fwd_pos > dst_rev_pos:
        mh_seq = str(fasta[dst_chr][dst_rev_pos:dst_fwd_pos]).upper()
        mh_len = dst_fwd_pos-dst_rev_pos
        dst_fwd_pos = dst_rev_pos

    src_chr, src_start, src_end = fwd_call.src_chr, fwd_call.src_pos, rev_call.src_pos
    if src_end-src_start < 50: continue

    seq = str(fasta[src_chr][src_start:src_end]).upper().replace("R", "N").replace("W", "N")
    if fwd_call.is_rc():
        seq = reverse_complement(seq)
    seq = mh_seq + seq

    precision = "IMPRECISE"
    if fwd_call.precision == "PRECISE" and rev_call.precision == "PRECISE":
        precision = "PRECISE"

    info = "MH_LEN=%d" % (mh_len)

    if len(seq) < 50 or dst_rev_pos-dst_fwd_pos >= len(seq): continue
    print(id, dst_chr, dst_fwd_pos, "N", dst_chr, dst_rev_pos, "N", "INS", seq, precision, info)
