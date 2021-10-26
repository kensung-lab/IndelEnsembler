from __future__ import print_function
import sys

def dir_to_strand(dir):
	if dir == 'F': return '+'
	else: return '-'

for line in sys.stdin:
	id, bp1, bp2 = line.split()[:3]
	dir1,chr1,pos1 = bp1[4:].split(':')
	dir2,chr2,pos2 = bp2[4:].split(':')
	id = id.replace("ID=", "SURVEYOR_")
	print(id, chr1, max(0,int(pos1)), dir_to_strand(dir1), chr2, max(0,int(pos2)), dir_to_strand(dir2), "TRA")

