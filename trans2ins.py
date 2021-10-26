from __future__ import print_function
import sys, pyfaidx
from collections import defaultdict

results = sys.stdin.readlines()
fasta = pyfaidx.Fasta(sys.argv[1])

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
def reverse_complement(seq):
	bases = list(seq)
	bases = reversed([complement.get(base,base) for base in bases])
	bases = ''.join(bases)
	return bases

ids = set()
both_sides_ids = set()
for line in results:
	id = line.split()[0]
	if id in ids:
		both_sides_ids.add(id)
	else:
		ids.add(id)

paired_calls = defaultdict(list)
for line in results:
	sl = line.split()
	if sl[0] in both_sides_ids:
		paired_calls[sl[0]].append(sl[1:])

for id in paired_calls:
	call1, call2 = paired_calls[id]

	dst_chr, dst_start, dst_end = call1[0], int(call1[1]), int(call2[1])
	if dst_start > dst_end: dst_start, dst_end = dst_end, dst_start

	src_chr, src_start, src_end = call1[3], int(call1[4]), int(call2[4])
	if src_start > src_end: src_start, src_end = src_end, src_start
	elif src_start == src_end: continue
	rc = call1[2] == call1[5]
	seq = str(fasta[src_chr][src_start:src_end]).upper()
	if rc:
		seq = reverse_complement(seq)

	print(id, dst_chr, dst_start, "N", dst_chr, dst_end, "N", "INS", seq)
