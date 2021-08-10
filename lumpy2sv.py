from __future__ import print_function
import sys, re, argparse, vcf

parser = argparse.ArgumentParser(description='Convert Lumpy output (stdin) into a standardized SV file.')
args = parser.parse_args()

vcf_reader = vcf.Reader(sys.stdin)
id = 0
for record in vcf_reader:
	svtype = record.INFO["SVTYPE"]

	if "SECONDARY" in record.INFO:
		continue

	if svtype == "BND":
		chr2_str = str(record.ALT[0])[1:-1].replace("[", "").replace("]", "").replace(":", " ")
		chr2, pos2 = chr2_str.split()
		pos2 = int(pos2)
		if record.CHROM == chr2:
			svtype = "INV"
		else:
			svtype = "TRA"
	else:
		chr2 = record.CHROM
		pos2 = int(record.INFO["END"])

	if svtype == "DEL":
		str1, str2 = "+", "-"
	elif svtype == "DUP":
		str1, str2 = "-", "+"
	elif svtype == "INV":
		str1, str2 = "N", "N"
	elif svtype == "TRA":
		strands = record.INFO["STRANDS"][0].split(":")[0]
		str1, str2 = strands

	if record.CHROM == chr2 and pos2-record.POS < 50: continue
	print("LUMPY_%d" %id, record.CHROM, record.POS, str1, chr2, pos2, str2, svtype)
	id += 1
