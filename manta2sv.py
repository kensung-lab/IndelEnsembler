from __future__ import print_function
import vcf, sys

vcf_reader = vcf.Reader(sys.stdin)
id = 0
for record in vcf_reader:
	if record.FILTER: continue

	svtype = record.INFO["SVTYPE"]
	if svtype == "BND":
		continue

	ins = ""
	if svtype == "INS" and ">" not in str(record.ALT[0]):
		ins = str(record.ALT[0])

	print("MANTA_%s"%id, record.CHROM, record.POS, "N", record.CHROM, record.INFO["END"], "N", svtype , ins)
	id += 1
