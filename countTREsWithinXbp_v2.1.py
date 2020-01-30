#!/usr/bin/python3

import os
import sys
import re
import argparse
import textwrap
import errno
from collections import defaultdict

def check_file(file):
	"""
	check that file exists
	"""

	if not os.path.isfile(file):
		print("Error! File does not exist:", file)
		print("Exiting now...")
		return False

	return True


def count_TREs(bed, gtf, bp, genebody, norm):
	"""
	count number of TREs within X bp from each gene/TSS
	"""

	# check that files exist
	bed_exist = check_file(bed)
	gtf_exist = check_file(gtf)
	if not (bed_exist and gtf_exist):
		return

	# check that files have the right extension
	if not bed.endswith(".bed"):
		print("Error! Bed file extension doesn't appear right:", bed)
		print("Exiting now...")
		return
	if not gtf.endswith(".gtf"):
		print("Error! GTF file extension doesn't appear right:", gtf)
		print("Exiting now...")
		return

	# read in gtf and store gene coordinates, strand information, and gene info in dictionary
	genes = defaultdict(lambda: defaultdict(dict))

	try:
		gtf_file = open(gtf)
	except:
		print("Error! Cannot open gtf file:", gtf)
		return

	for line in gtf_file:
		if line.strip().startswith("#"): # skip comments
			continue

		chrom, source, feature, start, end, score, strand, frame, attribute = line.split('\t')
		
		if not feature == "gene":
			continue

		# make sure chromosome formatting is "chr1" and not "1"
		chrom = str(chrom)
		if chrom.isdigit() or chrom in ["X", "Y", "MT", "M"]:
			chrom = "chr" + chrom

		string = attribute.split('"; ')
		ensemblpattern = re.compile("gene_id")
		ensembllist = list(filter(ensemblpattern.search, string))
		if len(ensembllist) != 1:
			print("Error! Maybe gtf attributes, especially gene id, are not formatted as expected!")
			print(ensembllist)
			return
		ensembl = ensembllist[0].split('id "')[1].split('"')[0]

		typepattern = re.compile("gene_.*type")
		typelist = list(filter(typepattern.search, string))
		if len(typelist) != 1:
			print("Error! Maybe gtf attributes, especially gene type, are not formatted as expected!")
			print(typelist)
			return
		genetype = typelist[0].split('type "')[1].split('"')[0]

		namepattern = re.compile("gene_name")
		namelist = list(filter(namepattern.search, string))
		if len(namelist) != 1:
			print("Error! Maybe gtf attributes, especially gene name, are not formatted as expected!")
			print(namelist)
			return
		gene = namelist[0].split('name "')[1].split('"')[0]
		gene = gene + "|" + ensembl + "|" + genetype
	
		if genebody: # window is around gene
			genes[chrom][gene]['lower'] = int(start) - bp - 1 # subtract 1 b/c bed is 0-based and gtf is 1-based coordinates
			genes[chrom][gene]['upper'] = int(end) + bp # don't subtract one here b/c bed is half-open and gtf is fully closed
		else: # window is around TSS
			if strand == "+":
				tss = int(start) - 1
			elif strand == "-":
				tss = int(end)
			else:
				print("Error! This gene strand is not + or -:", gene)
				return
			genes[chrom][gene]['lower'] = int(tss) - bp
			genes[chrom][gene]['upper'] = int(tss) + bp

		if genes[chrom][gene]['lower'] < 0:
			genes[chrom][gene]['lower'] = 0

		genes[chrom][gene]['count'] = 0

	gtf_file.close()
		
	# read in bed file and count TREs within X bp of gene ends

	try:
		bed_file = open(bed)
	except:
		print("Error! Cannot open bed file:", bed)
		return

	for line in bed_file:
		chromosome, first, last = line.split('\t')[:3]

		for key, value in genes[chromosome].items(): 
			if int(last) < genes[chromosome][key]['lower']:
				continue
			elif int(first) > genes[chromosome][key]['upper']:
				continue
			else:
				genes[chromosome][key]['count'] += 1

	bed_file.close()

	# print out gene, gene type, and # TREs within X bp
	for key, value in genes.items():
		for key2 in value.keys():
			if norm:
				outcount = (genes[key][key2]['count'] * 10000) / (genes[key][key2]['upper'] - genes[key][key2]['lower'])
			else:
				outcount = genes[key][key2]['count']
			try:
				sys.stdout.write(key2 + '\t' + str(outcount) + '\n')
			except IOError as e: #handle python pipe error
				if e.errno == errno.EPIPE:
					pass

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Count TREs Within X bp of each TSS/gene. Output is printed to stdout so users can pipe to additional functions and/or redirect to a file.", formatter_class=argparse.RawDescriptionHelpFormatter, epilog=textwrap.dedent('''Common pipes and redirects: 
Sort by number of TREs from high to low			| sort -k3,3nr
Extract all lincRNA (or other class) genes		| awk '$2 == "lincRNA" {print $0}'
Look at first 20 genes					| head -n 20
Redirect to a file					> filename.txt'''))
	parser.add_argument("bed", type=str, help="bed file of TREs")
	parser.add_argument("gtf", type=str, help="gencode gtf file")
	parser.add_argument("bp", type=int, help="number of bp from each TSS/gene to search for TREs. Must be integer")
	parser.add_argument("-g", "--gene", action="store_true", help="count TREs around gene instead of around TSS (Default: TSS)")
	parser.add_argument("-n", "--normalize", action="store_true", help="normalize number of TREs by genome search space. Useful for gene option because genes are different lengths. Automatically normalizes to TREs/10,000 bp")
	args = parser.parse_args()
	
	count_TREs(args.bed, args.gtf, args.bp, args.gene, args.normalize)
