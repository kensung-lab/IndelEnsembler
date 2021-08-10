import sys

def dir_to_strand(dir):
	if dir == 'R': return '+'
	if dir == 'L': return '-'
	return 'N'

for line in sys.stdin:
	id, bp1, bp2, svtype = line.split()[:4]
	size_str = line.split()[9].split('=')[1]
	if ':' in size_str:
		min_size, max_size = size_str.split(':')
		size = (int(min_size)+int(max_size))/2
	else:
		size = int(size_str)
	dir1,chr1,pos1 = bp1[4:].split(':')
	dir2,chr2,pos2 = bp2[4:].split(':')
	id = id.replace("ID=", "SURVEYOR_")
	print id, chr1, pos1, dir_to_strand(dir1), chr2, pos2, dir_to_strand(dir2), svtype, size
