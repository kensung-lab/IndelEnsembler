from __future__ import print_function
import sys

def dir_to_strand(dir):
    if dir == 'F': return '+'
    else: return '-'

for line in sys.stdin:
    sl = line.split()
    id, bp1, bp2 = sl[:3]
    sr = int(sl[4].split('=')[1])
    bp1_split, bp2_split = bp1[4:].split(':'), bp2[4:].split(':')
    dir1,chr1,pos1 = bp1_split[0], ":".join(bp1_split[1:-1]), bp1_split[-1]
    dir2,chr2,pos2 = bp2_split[0], ":".join(bp2_split[1:-1]), bp2_split[-1]

    precision = "IMPRECISE"
    if sr > 0:
        precision = "PRECISE"

    id = id.replace("ID=", "SURVEYOR_")
    print(id, chr1, max(0,int(pos1)), dir_to_strand(dir1), chr2, max(0,int(pos2)), dir_to_strand(dir2), "TRA", precision)
