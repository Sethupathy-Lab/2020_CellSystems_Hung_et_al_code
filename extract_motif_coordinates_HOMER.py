#!/usr/bin/python3

import os
import sys
import argparse
import textwrap
import errno
import re

def extract_motif(homer, chromNum, startNum, motifNum):
	"""
	extract motif sequences from HOMER annotatePeak.pl output
	"""
	
	for line in homer:
		line2 = line.strip().split('\t')
		chrom = line2[chromNum - 1]
		start = line2[startNum - 1]

		# skip header line if still there
		if chrom == "Chr":
			continue

		# skip lines without motifs if there are any or if not removed from file
		try:
			motif = line2[motifNum - 1]
		except IndexError as error:
			motif = None			
		if motif is None:
			continue

		# split motifs if there are multiple
		motif_list = motif.split("),")
		for single_motif in motif_list:
			motif_pos = single_motif.split("(")[0]
			motif_len = len(single_motif.split("(")[1].split(",")[0])
			motif_name = single_motif
			if motif_name.endswith(")"):
				pass
			else:
				motif_name = motif_name + ")"
			motif_name = str(chrom) + ":" + str(start) + ";" + motif_name
			strand = single_motif.split(",")[1]

			# print out motif coordinates for each motif
			startOut = int(start) + int(motif_pos) + 1
			stopOut = int(start) + int(motif_pos) + int(motif_len) + 1
			
			out = str(chrom) + "\t" + str(startOut) + "\t" + str(stopOut) + "\t" + motif_name + "\t" + str(0) + "\t" + str(strand)

			try:
				sys.stdout.write(out + '\n')
			except IOError as e: #handle python pipe error
				if e.errno == errno.EPIPE:
					pass

	homer.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Create bed file with motif coordinates from HOMER annotatePeak.pl output. Output is printed to stdout so users can pipe to additional functions and/or redirect to a file.") 
	parser.add_argument("homer", type=argparse.FileType('r'), help="homer annotatePeak.pl output. Must include at least peak coordinates and motif location within coordinates")
	parser.add_argument("-c", "--chrom", type=int, help="column of chr coordinate. 1-based (1st column is 1) (Default: 2)", default=2)
	parser.add_argument("-s", "--start", type=int, help="column of start coordinate (Default: 3)", default=3)
	parser.add_argument("-m", "--motif", type=int, help="column of motif locations (Default: 22)", default=22)
	args = parser.parse_args()
	
	extract_motif(args.homer, args.chrom, args.start, args.motif)
