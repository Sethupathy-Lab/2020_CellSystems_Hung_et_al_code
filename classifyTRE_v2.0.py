#!/usr/bin/python3

import os
import sys
import re
import argparse
import textwrap
import errno
from collections import defaultdict

def get_overlap(a, b):
	"""
	get overlap in bp of two intervals, a and b
	>0: number of bp overlap
	==0: adjacent to each other
	<0: distance apart
	"""
	return min(int(a[1]), int(b[1])) - max(int(a[0]), int(b[0]))


def classify_TRE(bed, gtf, down, up, overlap):
	"""
	count number of TREs within X bp from each gene
	"""
	
	# read in gtf and store gene coordinates and gene info in a dictionary and transcript coordinates and info in another dictionary 
	genes = defaultdict(lambda: defaultdict(dict))
	transcripts = defaultdict(lambda:defaultdict(dict))

	# loop thru gtf file and store TSS and gene features
	for line in gtf:
		if line.strip().startswith("#"): # skip comments
			continue

		chrom, source, feature, start, end, score, strand, frame, attribute = line.split('\t')
		
		if feature == "gene":

			# make sure chromosome formatting is "chr1" and not "1"
			chrom = str(chrom)
			if chrom.isdigit() or chrom in ["X", "Y", "MT", "M"]:
				chrom = "chr" + chrom

			string1 = attribute.split('"; ')
			typepattern = re.compile("gene_.*type")
			typelist = list(filter(typepattern.search, string1))
			if len(typelist) != 1:
				print("Error! Maybe gtf attributes, especially gene type, are not formatted as expected!")
				print(typelist)
				return
			genetype = typelist[0].split('type "')[1].split('"')[0]

			namepattern = re.compile("gene_name")
			namelist = list(filter(namepattern.search, string1))
			if len(namelist) != 1:
				print("Error! Maybe gtf attributes, especially gene name, are not formatted as expected!")
				print(namelist)
				return
			gene = namelist[0].split('name "')[1].split('"')[0]
	
			genes[chrom][gene]['lower'] = int(start) - 1 # subtract 1 b/c bed is 0-based and gtf is 1-based
			genes[chrom][gene]['upper'] = int(end) # don't subtract 1 b/c bed is half open and gtf is fully closed
			genes[chrom][gene]['type'] = genetype

		elif feature == "transcript":
			string1 = attribute.split('"; ')
			transcriptpattern = re.compile("transcript_id")
			transcriptlist = list(filter(transcriptpattern.search, string1))
			if len(transcriptlist) != 1:
				print("Error! Maybe gtf attributes, especially transcript id, are not formatted as expected!")
				print(transcriptlist)
				return
			transcript = transcriptlist[0].split('id "')[1].split('"')[0]

			typepattern = re.compile("gene_.*type")
			typelist = list(filter(typepattern.search, string1))
			if len(typelist) != 1:
				print("Error! Maybe gtf attributes, especially gene type, are not formatted as expected!")
				print(typelist)
				return
			genetype = typelist[0].split('type "')[1].split('"')[0]

			namepattern = re.compile("gene_name")
			namelist = list(filter(namepattern.search, string1))
			if len(namelist) != 1:
				print("Error! Maybe gtf attributes, especially gene name, are not formatted as expected!")
				print(namelist)
				return
			gene = namelist[0].split('name "')[1].split('"')[0]
			
			if strand == "+":
				transcripts[chrom][transcript]['lower'] = int(start) - 1 - up # transcript on plus strand so TSS is start
				transcripts[chrom][transcript]['upper'] = int(start) - 1 + down
			elif strand == "-":
				transcripts[chrom][transcript]['lower'] = int(end) - down # transcript on minus strand so TSS is end
				transcripts[chrom][transcript]['upper'] = int(end) + up 
			else:
				print("This transcript is not on any strand:", transcript)
				print("Skipping")
				continue
			transcripts[chrom][transcript]['type'] = genetype
			transcripts[chrom][transcript]['gene'] = gene

	gtf.close()

	# loop thru in bed file and classify each TRE
	for  line in bed:
		chromosome, first, last = line.split('\t')[:3]

		genestatus = False
		genelist = ""
		for key, value in genes[chromosome].items():
			over = get_overlap([genes[chromosome][key]['lower'], genes[chromosome][key]['upper']], [first, last])
			
			if over >= overlap: # bed region overlaps gene
				genelist = genelist + key + "," + genes[chromosome][key]['type'] + ";"
				genestatus = True
				
		transcriptstatus = False
		transcriptlist = ""
		for key, value in transcripts[chromosome].items():
			over = get_overlap([transcripts[chromosome][key]['lower'], transcripts[chromosome][key]['upper']], [first, last])

			if over >= overlap: # bed region overlaps promoter
				transcriptlist = transcriptlist + transcripts[chromosome][key]['gene'] + "," + key + ";"
				transcriptstatus = True

		if genestatus == False and transcriptstatus == False:
			distal = True
		else:
			distal = False

		if genelist == "":
			genelist = 'NA'

		if transcriptlist == "":
			transcriptlist = 'NA'

		geneprint = "\t" + "Gene:" + str(genestatus) + "\t" + genelist.rstrip(';')
		tssprint = "\t" + "Proximal:" + str(transcriptstatus) + "\t" + transcriptlist.rstrip(';')
		distalprint = "\t" + "Distal:" + str(distal) + "\n"

		try:
			sys.stdout.write(line.rstrip() + geneprint + tssprint + distalprint)
		except IOError as e: # handle python pipe error
			if e.errno == errno.EPIPE:
				pass

	bed.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Classify each TRE as distal (enhancer), proximal (promoter), or genic. TREs may be part of multiple genes (e.g. overlapping genes), multiple promoters (e.g. bidirectional genes), or both promoter and genic (e.g. multiple TSSs). All possibilities are output here. Output includes the TRE bed file followed by genic status, overlapping genes, proximal status, overlapping transcript TSSs, and distal status. Output is printed to stdout so users can pipe to additional functions and/or redirect to a file.", formatter_class=argparse.RawDescriptionHelpFormatter, epilog=textwrap.dedent('''Common pipes and redirects: 
Extract all proximal TREs				| awk '$8 == "Proximal:True" {print $0}'
Look at first 20 genes					| head -n 20
Redirect to a file					> filename.txt'''))
	parser.add_argument("bed", type=argparse.FileType('r'), help="bed file of TREs. Accepts a file or stdin", nargs="?", default = sys.stdin)
	parser.add_argument("gtf", type=argparse.FileType('r'), help="gencode gtf file")
	parser.add_argument("-u", "--up", type=int, help="distance in bp upstream of TSS to consider promoter (Default: 1000 bp)", default=1000)
	parser.add_argument("-d", "--down", type=int, help="distance in bp downstream of TSS to consider promoter (Default: 0 bp)", default=0)
	parser.add_argument("-o", "--overlap", type=int, help="bp of overlap required to be considered within a gene or promoter (Default: 1 bp)", default=1)	
	args = parser.parse_args()

	classify_TRE(args.bed, args.gtf, args.down, args.up, args.overlap)
